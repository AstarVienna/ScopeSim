# -*- coding: utf-8 -*-
"""
Effect for mapping spectral cubes to the detector plane.

The Effect is called `SpectralTraceList`, it applies a list of
`spectral_trace_list_utils.SpectralTrace` objects to a `FieldOfView`.
"""

from itertools import cycle
from typing import ClassVar

from tqdm.auto import tqdm

from astropy.io import fits
from astropy.table import Table

from .effects import Effect
from .ter_curves import FilterCurve
from .spectral_trace_list_utils import SpectralTrace, make_image_interpolations
from ..optics.image_plane_utils import header_from_list_of_xy
from ..base_classes import FieldOfViewBase, FOVSetupBase
from ..utils import from_currsys, check_keys, figure_factory, get_logger


logger = get_logger(__name__)


class SpectralTraceList(Effect):
    """
    List of spectral trace geometries for the detector plane.

    Should work in concert with an ApertureList (or ApertureMask) object and a
    DetectorList object

    Spectral trace patterns are to be kept in a ``fits.HDUList`` with one or
    more ``fits.BinTableHDU`` extensions, each one describing the geometry of a
    single trace. The first extension should be a ``BinTableHDU`` connecting
    the traces to the correct ``Aperture`` and ``ImagePlane`` objects.

    The ``fits.HDUList`` objects can be loaded using one of these two keywords:

    - ``filename``: for on disk FITS files, or
    - ``hdulist``: for in-memory ``fits.HDUList`` objects

    The format and contents of the extensions in the HDUList (FITS file) object
    is listed below

    Input Data Format
    -----------------
    A trace list FITS file needs the following extensions:

    - 0 : PrimaryHDU [header]
    - 1 : BinTableHDU [header, data] : Overview table of all traces
    - 2..N : BinTableHDU [header, data] : Trace tables. One per spectral trace

    EXT 0 : PrimaryHDU
    ++++++++++++++++++
    Required Header Keywords:

    - ECAT : int : Extension number of overview table. Normally 1
    - EDATA : int : Extension number of first Trace table. Normally 2

    No data is required in this extension

    EXT 1 : BinTableHDU : Overview of traces
    ++++++++++++++++++++++++++++++++++++++++
    No special header keywords are required in this extension

    Required Table columns:

    - description : str : identifier of each trace
    - extension_id : int : which extension is each trace in
    - aperture_id : int : which aperture matches this trace (e.g. MOS / IFU)
    - image_plane_id : int : on which image plane is this trace projected

    EXT 2 : BinTableHDU : Individual traces
    +++++++++++++++++++++++++++++++++++++++
    Required header keywords:
    - EXTNAME : must be identical to the `description` in EXT 1

    Recommended header keywords:
    - DISPDIR : "x" or "y" : dispersion axis. If not present, Scopesim tries
      to determine this automatically; this may be unreliable in some cases.

    Required Table columns:
    - wavelength : float : [um] : wavelength of monochromatic aperture image
    - s : float : [arcsec] : position along aperture perpendicular to trace
    - x : float : [mm] : x position of aperture image on focal plane
    - y : float : [mm] : y position of aperture image on focal plane

    """

    _class_params = {
        "x_colname": "x",
        "y_colname": "y",
        "s_colname": "s",
        "wave_colname": "wavelength",
        "col_number_start": 0,
        "center_on_wave_mid": False,
        "dwave": 0.002,  # [um] for finding best fit dispersion
        "invalid_value": None,  # for dodgy trace file values
    }
    z_order: ClassVar[tuple[int, ...]] = (70, 270, 670)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        params = {
            "pixel_scale": "!INST.pixel_scale",  # [arcsec / pix]}
            "plate_scale": "!INST.plate_scale",  # [arcsec / mm]
            "spectral_bin_width": "!SIM.spectral.spectral_bin_width", # [um]
            "wave_min": "!SIM.spectral.wave_min",  # [um]
            "wave_mid": "!SIM.spectral.wave_mid",  # [um]
            "wave_max": "!SIM.spectral.wave_max",  # [um]
            "x_colname": "x",
            "y_colname": "y",
            "s_colname": "s",
            "wave_colname": "wavelength",
            "center_on_wave_mid": False,
            "dwave": 0.002,  # [um] for finding the best fit dispersion
            "invalid_value": None,  # for dodgy trace file values
        }
        self.meta.update(params)

        # Parameters that are specific to the subclass
        self.meta.update(self._class_params)
        self.meta.update(kwargs)

        if self._file is not None:
            self.make_spectral_traces()

            self.update_meta()

    def make_spectral_traces(self):
        """Return a dictionary of spectral traces read in from a file."""
        self.ext_data = self._file[0].header["EDATA"]
        self.ext_cat = self._file[0].header["ECAT"]
        self.catalog = Table(self._file[self.ext_cat].data)
        spec_traces = {}
        for row in self.catalog:
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            hdu = self._file[row["extension_id"]]
            spec_traces[row["description"]] = SpectralTrace(hdu, **params)

        self.spectral_traces = spec_traces

    def update_meta(self):
        """
        Update fov related meta values.

        The values describe the full extent of the spectral trace
        volume in wavelength and space
        """
        wlim, xlim, ylim = [], [], []
        for thetrace in self.spectral_traces.values():
            fov = thetrace.fov_grid()
            if "wave_min" in fov:
                wlim.extend([fov["wave_min"], fov["wave_max"]])
            if "x_min" in fov:
                xlim.extend([fov["x_min"], fov["x_max"]])
            if "y_min" in fov:
                ylim.extend([fov["y_min"], fov["y_max"]])

        if wlim:
            self.meta["wave_min"] = min(wlim)
            self.meta["wave_max"] = max(wlim)
        if xlim:
            self.meta["x_min"] = min(xlim)
            self.meta["x_max"] = max(xlim)
        if ylim:
            self.meta["y_min"] = min(ylim)
            self.meta["y_max"] = max(ylim)

    def apply_to(self, obj, **kwargs):
        """
        Interface between ``FieldOfView`` and ``SpectralTraceList``.

        This is called twice:
        1. During setup of the required FieldOfView objects, the
        SpectralTraceList is asked for the source space volumes that
        it requires (spatial limits and wavelength limits).
        2. During "observation" the method is passed a single FieldOfView
        object and applies the mapping to the image plane to it.
        The FieldOfView object is associated to one SpectralTrace from the
        list, identified by meta["trace_id"].
        """
        if isinstance(obj, FOVSetupBase):
            # Setup of FieldOfView object
            # volumes = [spectral_trace.fov_grid()
            #            for spectral_trace in self.spectral_traces.values()]

            new_vols_list = []

            # for vol in volumes:
            for spt in self.spectral_traces.values():
                vol = spt.fov_grid()
                wave_edges = [vol["wave_min"], vol["wave_max"]]
                if "x_min" in vol:
                    x_edges = [vol["x_min"], vol["x_max"]]
                    y_edges = [vol["y_min"], vol["y_max"]]
                    extracted_vols = obj.extract(
                        axes=["wave", "x", "y"],
                        edges=(wave_edges, x_edges, y_edges),
                        aperture_id=vol["aperture_id"])
                else:
                    extracted_vols = obj.extract(
                        axes=["wave"],
                        edges=(wave_edges, ),
                        aperture_id=vol["aperture_id"])

                for ex_vol in extracted_vols:
                    ex_vol["meta"].update(vol)
                    ex_vol["meta"].pop("wave_min")
                    ex_vol["meta"].pop("wave_max")
                new_vols_list.extend(extracted_vols)

            obj.volumes = new_vols_list

        if isinstance(obj, FieldOfViewBase):
            # Application to field of view
            if obj.hdu is not None and obj.hdu.header["NAXIS"] == 3:
                obj.cube = obj.hdu
            elif obj.hdu is not None and obj.hdu.header["NAXIS"] == 2:
                # todo: catch the case of obj.hdu.header["NAXIS"] == 2
                # for MAAT
                pass
            elif obj.hdu is None and obj.cube is None:
                logger.info("Making cube")
                obj.cube = obj.make_cube_hdu()

            spt = self.spectral_traces[obj.trace_id]
            obj.hdu = spt.map_spectra_to_focal_plane(obj)

        return obj

    @property
    def footprint(self):
        """Return the footprint of the entire SpectralTraceList."""
        xfoot, yfoot = [], []
        for spt in self.spectral_traces.values():
            xtrace, ytrace = spt.footprint()
            xfoot.extend(xtrace)
            yfoot.extend(ytrace)

        xfoot = [min(xfoot), max(xfoot), max(xfoot), min(xfoot)]
        yfoot = [min(yfoot), min(yfoot), max(yfoot), max(yfoot)]

        return xfoot, yfoot

    @property
    def image_plane_header(self):
        """Create and return header for the ImagePlane."""
        x, y = self.footprint
        pixel_scale = from_currsys(self.meta["pixel_scale"], self.cmds)
        hdr = header_from_list_of_xy(x, y, pixel_scale, "D")

        return hdr

    def rectify_traces(self, hdulist, xi_min=None, xi_max=None, interps=None,
                       **kwargs):
        """Create rectified 2D spectra for all traces in the list.

        This method creates an HDU list with one extension per spectral
        trace, i.e. it essentially treats all traces independently.
        For the case of an IFU where the traces correspond to spatial
        slices for the same wavelength range, use method `rectify_cube`
        (not yet implemented).

        Parameters
        ----------
        hdulist : str or fits.HDUList
           The result of scopesim readout()
        xi_min, xi_max : float [arcsec]
           Spatial limits of the slit on the sky. This should be taken
           from the header of the hdulist, but this is not yet provided by
           scopesim. For the time being, these limits *must* be provided by
           the user.
        interps :  list of interpolation functions
           If provided, there must be one for each image extension in
           `hdulist`. The functions go from pixels to the images and can be
           created with, e.g. ``RectBivariateSpline``.
        """
        try:
            inhdul = fits.open(hdulist)
        except TypeError:
            inhdul = hdulist

        # Crude attempt to get a useful wavelength range
        # Problematic because different instruments use different
        # keywords for the filter... We try to make it work for METIS
        # and MICADO for the time being.
        try:
            filter_name = from_currsys("!OBS.filter_name", self.cmds)
        except ValueError:
            filter_name = from_currsys("!OBS.filter_name_fw1", self.cmds)

        filtcurve = FilterCurve(
            filter_name=filter_name,
            filename_format=from_currsys("!INST.filter_file_format", self.cmds))
        filtwaves = filtcurve.table["wavelength"]
        filtwave = filtwaves[filtcurve.table["transmission"] > 0.01]
        wave_min, wave_max = min(filtwave), max(filtwave)
        logger.info(
            "Full wavelength range: %.02f .. %.02f um", wave_min, wave_max)

        if xi_min is None or xi_max is None:
            try:
                xi_min = inhdul[0].header["HIERARCH INS SLIT XIMIN"]
                xi_max = inhdul[0].header["HIERARCH INS SLIT XIMAX"]
                logger.info(
                    "Slit limits taken from header: %.02f .. %.02f arcsec",
                    xi_min, xi_max)
            except KeyError:
                logger.error(
                    "Spatial slit limits (in arcsec) must be provided:\n"
                    "- either as method parameters xi_min and xi_max\n"
                    "- or as header keywords HIERARCH INS SLIT XIMIN/XIMAX"
                )
                return None

        bin_width = kwargs.get("bin_width", None)

        if interps is None:
            logger.info("Computing interpolation functions")
            interps = make_image_interpolations(hdulist)

        pdu = fits.PrimaryHDU()
        pdu.header["FILETYPE"] = "Rectified spectra"
        # pdu.header["INSTRUME"] = inhdul[0].header["HIERARCH ESO OBS INSTRUME"]
        # pdu.header["FILTER"] = from_currsys("!OBS.filter_name_fw1", self.cmds)
        outhdul = fits.HDUList([pdu])

        for i, trace_id in tqdm(enumerate(self.spectral_traces, start=1),
                                desc=" Traces"):
            hdu = self[trace_id].rectify(hdulist,
                                         interps=interps,
                                         bin_width=bin_width,
                                         xi_min=xi_min, xi_max=xi_max,
                                         wave_min=wave_min, wave_max=wave_max)
            if hdu is not None:   # ..todo: rectify does not do that yet
                outhdul.append(hdu)
                outhdul[0].header[f"EXTNAME{i}"] = trace_id

        outhdul[0].header.update(inhdul[0].header)

        return outhdul

    def rectify_cube(self, hdulist):
        """Rectify traces and combine into a cube."""
        raise NotImplementedError()

    def plot(self, wave_min=None, wave_max=None, axes=None, **kwargs):
        """Plot every spectral trace in the spectral trace list.

        Parameters
        ----------
        wave_min : float, optional
            Minimum wavelength, if any. If None, value from_currsys is used.
        wave_max : float, optional
            Maximum wavelength, if any. If None, value from_currsys is used.
        axes : matplotlib axes, optional
            The axes object to use for the plot. If None (default), a new
            figure with one axes will be created.
        **kwargs : dict
            Any other parameters passed along to the plot method of the
            individual spectral traces.

        Returns
        -------
        fig : matplotlib figure
            DESCRIPTION.

        """
        if wave_min is None:
            wave_min = from_currsys("!SIM.spectral.wave_min", self.cmds)
        if wave_max is None:
            wave_max = from_currsys("!SIM.spectral.wave_max", self.cmds)

        if axes is None:
            fig, axes = figure_factory()
        else:
            fig = axes.figure

        if self.spectral_traces is not None:
            for spt, c in zip(self.spectral_traces.values(), cycle("rgbcymk")):
                spt.plot(wave_min, wave_max, c=c, axes=axes, **kwargs)

        return fig

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(**{self.meta!r})"

    def __str__(self) -> str:
        msg = (f"SpectralTraceList: \"{self.meta.get('name')}\" : "
               f"{len(self.spectral_traces)} traces")
        return msg

    def __getitem__(self, item):
        return self.spectral_traces[item]

    def __setitem__(self, key, value):
        self.spectral_traces[key] = value


class SpectralTraceListWheel(Effect):
    """
    A Wheel-Effect object for selecting between multiple gratings/grisms.

    See ``SpectralTraceList`` for the trace file format description.

    Parameters
    ----------
    trace_list_names : list
        The list of unique identifiers in the trace filenames

    filename_format : str
        ``f-string`` that directs scopesim to the folder containing the trace
        files. This can be a ``!-string`` if the trace names are shared with
        other ``*Wheel`` effect objects (e.g. a ``FilterWheel``). See examples.

    current_trace_list : str
        default trace file to use

    kwargs : key-value pairs
        Addition keywords that are passed to the ``SpectralTraceList`` objects
        See SpectralTraceList docstring

    Examples
    --------
    A simplified YAML file example taken from the OSIRIS instrument package::

        alias: INST
        name: OSIRIS_LSS

        properties:
          decouple_detector_from_sky_headers: True
          grism_names:
            - R300B
            - R500B
            - R1000B
            - R2500V

        effects:
          - name: spectral_trace_wheel
            description: grism wheel contining spectral trace geometries
            class: SpectralTraceListWheel
            kwargs:
              current_trace_list: "!OBS.grating_name"
              filename_format: "traces/LSS_{}_TRACE.fits"
              trace_list_names: "!INST.grism_names"

          - name: grating_efficiency
            description: OSIRIS grating efficiency curves, piggybacking on FilterWheel
            class: FilterWheel
            kwargs:
              minimum_throughput: !!float 0.
              filename_format: "gratings/{}.txt"
              current_filter: "!OBS.grating_name"
              filter_names: "!INST.grism_names"

    """

    required_keys = {
        "trace_list_names",
        "filename_format",
        "current_trace_list",
    }
    z_order: ClassVar[tuple[int, ...]] = (70, 270, 670)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = True
    report_table_rounding: ClassVar[int] = 4
    _current_str = "current_trace_list"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(kwargs, self.required_keys, action="error")

        params = {
            "path": "",
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        path = self._get_path()
        self.trace_lists = {}
        if "name" in kwargs:
            kwargs.pop("name")
        for name in from_currsys(self.meta["trace_list_names"], self.cmds):
            fname = str(path).format(name)
            self.trace_lists[name] = SpectralTraceList(filename=fname,
                                                       name=name,
                                                       **kwargs)

    def apply_to(self, obj, **kwargs):
        """Use apply_to of current trace list."""
        return self.current_trace_list.apply_to(obj, **kwargs)

    @property
    def current_trace_list(self):
        trace_list_eff = None
        trace_list_name = from_currsys(self.meta["current_trace_list"],
                                       self.cmds)
        if trace_list_name is not None:
            trace_list_eff = self.trace_lists[trace_list_name]
        return trace_list_eff
