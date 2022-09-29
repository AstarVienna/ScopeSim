"""
Effect for mapping spectral cubes to the detector plane

The Effect is called SpectralTraceList, it applies a list of
optics.spectral_trace_SpectralTrace objects to a FieldOfView.
"""

from os import path as pth
import numpy as np

from astropy.io import fits
from astropy.table import Table

from .effects import Effect
from .spectral_trace_list_utils import SpectralTrace
from ..utils import from_currsys, check_keys, interp2
from ..optics.image_plane_utils import header_from_list_of_xy
from ..base_classes import FieldOfViewBase, FOVSetupBase


class SpectralTraceList(Effect):
    """
    List of spectral trace geometries for the detector plane

    Should work in concert with an ApertureList (or ApertureMask) object and a
    DetectorList object

    Spectral trace patterns are to be kept in a ``fits.HDUList`` with one or
    more ``fits.BinTableHDU`` extensions, each one describing the geometry of a
    single trace. The first extension should be a ``BinTableHDU`` connecting the
    traces to the correct ``Aperture`` and ``ImagePlane`` objects.

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

    - description : str : description of each each trace
    - extension_id : int : which extension is each trace in
    - aperture_id : int : which aperture matches this trace (e.g. MOS / IFU)
    - image_plane_id : int : on which image plane is this trace projected

    EXT 2 : BinTableHDU : Individual traces
    +++++++++++++++++++++++++++++++++++++++
    No special header keywords are required in this extension

    Required Table columns:

    - wavelength : float : [um] : wavelength of monochromatic aperture image
    - s : float : [arcsec] : position along aperture perpendicular to trace
    - x : float : [mm] : x position of aperture image on focal plane
    - y : float : [mm] : y position of aperture image on focal plane

    """
    _class_params = {"x_colname": "x",
                     "y_colname": "y",
                     "s_colname": "s",
                     "wave_colname": "wavelength",
                     "col_number_start": 0,
                     "center_on_wave_mid": False,
                     "dwave": 0.002,  # [um] for finding the best fit dispersion
                     "invalid_value": None,  # for dodgy trace file values
                     }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        params = {"z_order": [70, 270, 670],
                  "pixel_scale": "!INST.pixel_scale",  # [arcsec / pix]}
                  "plate_scale": "!INST.plate_scale",  # [arcsec / mm]
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
                  "report_plot_include": True,
                  "report_table_include": False,
                  }
        self.meta.update(params)

        # Parameters that are specific to the subclass
        self.meta.update(self._class_params)
        self.meta.update(kwargs)

        if self._file is not None:
            self.make_spectral_traces()

    def make_spectral_traces(self):
        """Returns a dictionary of spectral traces read in from a file"""
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

    def apply_to(self, obj, **kwargs):
        """
        Interface between FieldOfView and SpectralTraceList

        This is called twice:
        1. During setup of the required FieldOfView objects, the
        SpectralTraceList is asked for the source space volumes that
        it requires (spatial limits and wavelength limits).
        2. During "observation" the method is passed a single FieldOfView
        object and applies the mapping to the image plane to it.
        The FieldOfView object is associated to one SpectralTrace from the
        list, identified by meta['trace_id'].
        """
        if isinstance(obj, FOVSetupBase):
            # Setup of FieldOfView object
            volumes = [self.spectral_traces[key].fov_grid()
                       for key in self.spectral_traces]
            new_vols_list = []
            for vol in volumes:
                wave_edges = [vol["wave_min"], vol["wave_max"]]
                if "x_min" in vol:
                    x_edges = [vol["x_min"], vol["x_max"]]
                    y_edges = [vol["y_min"], vol["y_max"]]
                    extracted_vols = obj.extract(axes=["wave", "x", "y"],
                                                 edges=(wave_edges,
                                                        x_edges,
                                                        y_edges),
                                                 aperture_id=vol["aperture_id"])
                else:
                    extracted_vols = obj.extract(axes=["wave"],
                                                 edges=(wave_edges, ),
                                                 aperture_id=vol["aperture_id"])

                for ex_vol in extracted_vols:
                    ex_vol["meta"].update(vol)
                    ex_vol["meta"].pop("wave_min")
                    ex_vol["meta"].pop("wave_max")
                new_vols_list += extracted_vols

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
                obj.cube = obj.make_cube_hdu()

            # ..todo: obj will be changed to a single one covering the full field of view
            # covered by the image slicer (28 slices for LMS; for LSS still only a single slit)
            # We need a loop over spectral_traces that chops up obj into the single-slice fov before
            # calling map_spectra...
            trace_id = obj.meta['trace_id']
            spt = self.spectral_traces[trace_id]
            obj.hdu = spt.map_spectra_to_focal_plane(obj)

        return obj

    def get_waveset(self, pixel_size=None):
        if pixel_size is None:
            pixel_size = self.meta["pixel_scale"] / self.meta["plate_scale"]

        wavesets = [spt.get_pixel_wavelength_edges(pixel_size)
                    for spt in self.spectral_traces]

        return wavesets

    def get_fov_headers(self, sky_header, **kwargs):
        fov_headers = []
        for spt in self.spectral_traces:
            fov_headers += spt.fov_headers(sky_header=sky_header, **kwargs)

        return fov_headers


    @property
    def footprint(self):
        """Return the footprint of the entire SpectralTraceList"""
        xfoot, yfoot = [], []
        for spt in self.spectral_traces.values():
            xtrace, ytrace = spt.footprint()
            xfoot += xtrace
            yfoot += ytrace

        xfoot = [np.min(xfoot), np.max(xfoot), np.max(xfoot), np.min(xfoot)]
        yfoot = [np.min(yfoot), np.min(yfoot), np.max(yfoot), np.max(yfoot)]

        return xfoot, yfoot

    @property
    def image_plane_header(self):
        x, y = self.footprint
        pixel_scale = from_currsys(self.meta["pixel_scale"])
        hdr = header_from_list_of_xy(x, y, pixel_scale, "D")

        return hdr

    def plot(self, wave_min=None, wave_max=None, **kwargs):
        if wave_min is None:
            wave_min = from_currsys("!SIM.spectral.wave_min")
        if wave_max is None:
            wave_max = from_currsys("!SIM.spectral.wave_max")

        from matplotlib import pyplot as plt
        from matplotlib._pylab_helpers import Gcf
        if len(Gcf.figs) == 0:
            plt.figure(figsize=(12, 12))

        if self.spectral_traces is not None:
            clrs = "rgbcymk" * (1 + len(self.spectral_traces) // 7)
            for spt, c in zip(self.spectral_traces.values(), clrs):
                spt.plot(wave_min, wave_max, c=c)

        return plt.gcf()

    def __repr__(self):
        return "\n".join([spt.__repr__() for spt in self.spectral_traces])

    def __str__(self):
        msg = 'SpectralTraceList: "{}" : {} traces' \
              ''.format(self.meta.get("name"), len(self.spectral_traces))
        return msg


class SpectralTraceListWheel(Effect):
    """
    A Wheel-Effect object for selecting between multiple gratings/grisms

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
    def __init__(self, **kwargs):
        required_keys = ["trace_list_names",
                         "filename_format",
                         "current_trace_list"]
        check_keys(kwargs, required_keys, action="error")

        super(SpectralTraceListWheel, self).__init__(**kwargs)

        params = {"z_order": [70, 270, 670],
                  "path": "",
                  "report_plot_include": True,
                  "report_table_include": True,
                  "report_table_rounding": 4}
        self.meta.update(params)
        self.meta.update(kwargs)

        path = pth.join(self.meta["path"],
                        from_currsys(self.meta["filename_format"]))
        self.trace_lists = {}
        for name in from_currsys(self.meta["trace_list_names"]):
            kwargs["name"] = name
            self.trace_lists[name] = SpectralTraceList(filename=path.format(name),
                                                       **kwargs)

    def apply_to(self, obj, **kwargs):
        """Use apply_to of current trace list"""
        return self.current_trace_list.apply_to(obj, **kwargs)

    @property
    def current_trace_list(self):
        trace_list_eff = None
        trace_list_name = from_currsys(self.meta["current_trace_list"])
        if trace_list_name is not None:
            trace_list_eff = self.trace_lists[trace_list_name]
        return trace_list_eff

    @property
    def display_name(self):
        name = self.meta.get("name", self.meta.get("filename", "<untitled>"))
        return f'{name} : [{from_currsys(self.meta["current_trace_list"])}]'