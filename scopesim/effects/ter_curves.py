# -*- coding: utf-8 -*-
"""Transmission, emissivity, reflection curves."""

import warnings
from typing import ClassVar
from collections.abc import Collection, Iterable

import numpy as np
import skycalc_ipy
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

from .effects import Effect
from .ter_curves_utils import (add_edge_zeros, combine_two_spectra,
                               apply_throughput_to_cube, download_svo_filter,
                               download_svo_filter_list)
from ..base_classes import SourceBase, FOVSetupBase
from ..optics.surface import SpectralSurface
from ..source.source import Source
from ..source.source_fields import CubeSourceField, SpectrumSourceField
from ..utils import (from_currsys, quantify, check_keys, find_file,
                     figure_factory, get_logger)


logger = get_logger(__name__)


class TERCurve(Effect):
    """
    Transmission, Emissivity, Reflection Curve.

    note:: This is basically an ``Effect`` wrapper for the
           ``SpectralSurface`` object

    Must contain a wavelength column, and one or more of the following:
    ``transmission``, ``emissivity``, ``reflection``.
    Additionally, in the header there
    should be the following keywords: wavelength_unit

    kwargs that can be passed::

        "rescale_emission" : { "filter_name": str, "value": float, "unit": str}

    Examples
    --------
    Directly inside a YAML file description::

        name: bogus_surface
        class: TERCurve
        kwargs:
            array_dict:
                wavelength: [0.3, 3.0]
                transmission: [0.9, 0.9]
                emission: [1, 1]
            wavelength_unit: um
            emission_unit: ph s-1 m-2 um-1
            rescale_emission:
                filter_name: "Paranal/HAWK.Ks"
                value: 15.5
                unit: ABmag

    Indirectly inside a YAML file::

        name: some_curve
        class TERCurve
        kwargs:
            filename: bogus_surface.dat

    which references this ASCII file::

        # name: bogus_surface
        # wavelength_unit: um
        wavelength  transmission    emissivity
        0.3         0.9             0.1
        3.0         0.9             0.1

    """

    z_order: ClassVar[tuple[int, ...]] = (10, 110, 510)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, filename=None, **kwargs):
        super().__init__(filename=filename, **kwargs)
        params = {
            "ignore_wings": False,
            "wave_min": "!SIM.spectral.wave_min",
            "wave_max": "!SIM.spectral.wave_max",
            "wave_unit": "!SIM.spectral.wave_unit",
            "wave_bin": "!SIM.spectral.spectral_bin_width",
            "bg_cell_width": "!SIM.computing.bg_cell_width",
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        self.surface = SpectralSurface(cmds=self.cmds)
        self.surface.meta.update(self.meta)
        self._background_source = None

        data = self.data
        if self.meta["ignore_wings"]:
            data = add_edge_zeros(data, "wavelength")
        if data is not None:
            # Assert that get_data() did not give us an image.
            assert isinstance(data, Table), "TER Curves must be tables."
            self.surface.table = data
            self.surface.table.meta.update(self.meta)

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, SourceBase):
            assert isinstance(obj, Source), "Only Source supported."
            self.meta = from_currsys(self.meta, self.cmds)
            wave_min = quantify(self.meta["wave_min"], u.um).to(u.AA)
            wave_max = quantify(self.meta["wave_max"], u.um).to(u.AA)

            thru = self.throughput

            # apply transmission to source spectra
            for fld in obj.fields:
                if isinstance(fld, CubeSourceField):
                    fld.field = apply_throughput_to_cube(
                        fld.field, thru, fld.wave)
                elif isinstance(fld, SpectrumSourceField):
                    fld.spectra = {
                        isp: combine_two_spectra(spec, thru, "multiply",
                                                 wave_min, wave_max)
                        for isp, spec in fld.spectra.items()
                    }
                else:
                    # Rather log than raise here, can still move on
                    logger.error("Source field is neither Cube nor has "
                                 "spectra, this shouldn't occur...")

            # add the effect background to the source background field
            if self.background_source is not None:
                obj.append(self.background_source)

        if isinstance(obj, FOVSetupBase):
            from ..optics.fov_manager import FovVolumeList
            assert isinstance(obj, FovVolumeList), "Only FovVolumeList supported."
            wave = self.surface.throughput.waveset
            thru = self.surface.throughput(wave)
            valid_waves = np.argwhere(thru > 0)
            wave_min = wave[max(0, valid_waves[0][0] - 1)]
            wave_max = wave[min(len(wave) - 1, valid_waves[-1][0] + 1)]

            obj.shrink("wave", [wave_min.to(u.um).value,
                                wave_max.to(u.um).value])

        return obj

    @property
    def emission(self):
        return self.surface.emission

    @property
    def throughput(self):
        return self.surface.throughput

    @property
    def background_source(self):
        if self._background_source is None:
            flux = self.emission
            bg_hdu = fits.ImageHDU()

            bg_hdu.header.update({"BG_SRC": True,
                                  "BG_SURF": self.display_name,
                                  "CUNIT1": "ARCSEC",
                                  "CUNIT2": "ARCSEC",
                                  "CDELT1": 0,
                                  "CDELT2": 0,
                                  "BUNIT": "PHOTLAM arcsec-2",
                                  "SOLIDANG": "arcsec-2"})
            self._background_source = Source(image_hdu=bg_hdu, spectra=flux)

        return self._background_source

    def plot(self, which="x", wavelength=None, *, axes=None, **kwargs):
        """Plot TER curves.

        Parameters
        ----------
        which : {"x", "t", "e", "r"}, optional
            "x" plots throughput. "t","e","r" plot trans/emission/refl.
            Can be a combination, e.g. "tr" or "tex" to plot each.
        wavelength : array_like, optional
            Wavelength on x-axis, taken from currsys if None (default).
        axes : matplotlib axes, optional
            If given, plot into existing axes. The default is None.

        Returns
        -------
        fig : matplotlib figure
            Figure containing plots.

        """
        if axes is None:
            fig, axes = figure_factory(len(which), 1, iterable_axes=True)
        else:
            fig = axes.figure
            _guard_plot_axes(which, axes)

        self.meta.update(kwargs)
        params = from_currsys(self.meta, self.cmds)

        wave_unit = self.meta.get("wavelength_unit")
        if wavelength is None:
            wunit = params["wave_unit"]
            # TODO: shouldn't need both, make sure they're equal
            if wunit != wave_unit:
                logger.warning("wavelength units in the meta dict of "
                             "%s are inconsistent:\n"
                             "- wavelength_unit : %s\n"
                             "- wave_unit : %s",
                             {self.meta.get("name")},
                             wave_unit, wunit)

            wave = np.arange(quantify(params["wave_min"], wunit).value,
                             quantify(params["wave_max"], wunit).value,
                             quantify(params["wave_bin"], wunit).value)
            wave *= u.Unit(wunit)
        else:
            wave = wavelength

        plot_kwargs = self.meta.get("plot_kwargs", {})
        abbrs = {"t": "transmission", "e": "emission",
                 "r": "reflection", "x": "throughput"}

        if not isinstance(axes, Iterable):
            axes = [axes]
        for ter, ax in zip(which, axes):
            y_name = abbrs.get(ter, "throughput")
            y = getattr(self.surface, y_name)
            if not isinstance(y, u.Quantity):  # assume synphot spectrum
                y = y(wave)
            ax.plot(wave, y, **plot_kwargs)
            ax.set_xlabel(f"Wavelength [{wave_unit}]")
            y_unit = str(y.unit) or "dimensionless"
            ax.set_ylabel(f"{y_name.title()} [{y_unit}]")

        return fig


class AtmosphericTERCurve(TERCurve):
    z_order: ClassVar[tuple[int, ...]] = (111, 511)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["action"] = "transmission"
        self.meta["position"] = 0       # position in surface table
        self.meta.update(kwargs)
        self.surface.meta.update(self.meta)


class SkycalcTERCurve(AtmosphericTERCurve):
    """
    Retrieve an atmospheric spectrum from ESO's skycalc server.

    kwarg parameters
    ----------------
    skycalc parameters can be found by calling::

        >>> import skycalc_ipy
        >>> skycalc_ipy.SkyCalc().keys

    .. note:: Different to ``skycalc_ipy``, `wmin` and `wmax` must be given in
        units of ``um``

    Examples
    --------
    ::

        - name : skycalc_background
          class : SkycalcTERCurve
          kwargs :
            wunit : "!SIM.spectral.wave_unit"
            wmin : "!SIM.spectral.wave_min"
            wmax : "!SIM.spectral.wave_max"
            wdelta : 0.0001     # 0.1nm bin width
            outer : 1
            outer_unit : "m"

    """

    z_order: ClassVar[tuple[int, ...]] = (112, 512)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["use_local_skycalc_file"] = False
        self.meta.update(kwargs)

        self.skycalc_table = None
        self.skycalc_conn = None

        if self.include is True:
            # Only query the database if the effect is actually included.
            # Sets skycalc_conn and skycalc_table.
            self.load_skycalc_table()

    @property
    def include(self):
        return from_currsys(self.meta["include"], self.cmds)

    @include.setter
    def include(self, item):
        self.meta["include"] = item
        if item is True and self.skycalc_table is None:
            self.load_skycalc_table()

    def load_skycalc_table(self):
        use_local_file = from_currsys(self.meta["use_local_skycalc_file"],
                                      self.cmds)
        if not use_local_file:
            self.skycalc_conn = skycalc_ipy.SkyCalc()
            tbl = self.query_server()

            if "name" not in self.meta:
                self.meta["name"] = self.skycalc_conn["observatory"]

        else:
            path = find_file(use_local_file)
            fits_tbl = fits.getdata(path, ext=1)
            fits_hdr = fits.getheader(path, ext=0)
            tbl = Table(fits_tbl)
            tbl["lam"].unit = u.um
            for colname in tbl.colnames:
                if "flux" in colname:
                    tbl[colname].unit = u.Unit("ph s-1 m-2 um-1 arcsec-2")
            tbl_small = Table()
            tbl_small.meta["fits_header"] = dict(fits_hdr)
            tbl_small.add_columns([tbl["lam"], tbl["trans"], tbl["flux"]])
            tbl = tbl_small

        for i, colname in enumerate(["wavelength", "transmission", "emission"]):
            tbl.columns[i].name = colname
        tbl.meta["wavelength_unit"] = tbl.columns[0].unit
        tbl.meta["emission_unit"] = tbl.columns[2].unit

        self.surface.table = tbl
        self.surface.meta.update(tbl.meta)
        self.skycalc_table = tbl

    def query_server(self, **kwargs):
        self.meta.update(kwargs)

        if "wunit" in self.meta:
            scale_factor = u.Unit(from_currsys(self.meta["wunit"],
                                               self.cmds)).to(u.nm)
            for key in ["wmin", "wmax", "wdelta"]:
                if key in self.meta:
                    self.meta[key] = from_currsys(self.meta[key],
                                                  self.cmds) * scale_factor

        conn_kwargs = {key: self.meta[key] for key in self.meta
                       if key in self.skycalc_conn.defaults}
        conn_kwargs = from_currsys(conn_kwargs, self.cmds)
        self.skycalc_conn.values.update(conn_kwargs)

        try:
            tbl = self.skycalc_conn.get_sky_spectrum(return_type="table")
        except ConnectionError:
            msg = "Could not connect to skycalc server"
            logger.exception(msg)
            raise ValueError(msg)

        return tbl


class QuantumEfficiencyCurve(TERCurve):
    z_order: ClassVar[tuple[int, ...]] = (113, 513)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["action"] = "transmission"
        self.meta["position"] = -1          # position in surface table


class FilterCurve(TERCurve):
    """
    Descripton TBA.

    Parameters
    ----------
    position : int, optional
    filter_name : str, optional
        "Ks" - corresponding to the filter name in the filename pattern
    filename_format : str, optional
        "TC_filter_{}.dat"

    Can either be created using the standard 3 options:
    - `filename`: direct filename of the filter curve
    - `table`: an ``astropy.Table``
    - `array_dict`: a dictionary version of a table: ``{col_name1: values, }``

    or by passing the combination of `filter_name` and `filename_format` as
    kwargs. Here all filter file names follow a pattern (e.g. see above) and
    the "{}" are replaced by `filter_name` at run time. `filter_name` can
    also be a !-string for a ``__currsys__`` entry, e.g. "!INST.filter_name".

    """

    z_order: ClassVar[tuple[int, ...]] = (114, 214, 514)

    def __init__(self, cmds=None, **kwargs):
        # super().__init__(**kwargs)
        if not np.any([key in kwargs for key in ["filename", "table",
                                                 "array_dict"]]):
            if "filter_name" in kwargs and "filename_format" in kwargs:
                filt_name = from_currsys(kwargs["filter_name"], cmds)
                file_format = from_currsys(kwargs["filename_format"], cmds)
                kwargs["filename"] = file_format.format(filt_name)
            else:
                raise ValueError("FilterCurve must be passed one of "
                                 "(`filename`, `array_dict`, `table`) or both "
                                 f"(`filter_name`, `filename_format`): {kwargs}")

        super().__init__(cmds=cmds, **kwargs)
        if self.table is None:
            raise ValueError("Could not initialise filter. Either filename "
                             "not found, or array are not compatible")

        params = {"minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "action": "transmission",
                  "position": -1,               # position in surface table
                  "wing_flux_level": None,
                  "name": "untitled filter"}
        self.meta.update(params)
        self.meta.update(kwargs)

        min_thru = from_currsys(self.meta["minimum_throughput"], self.cmds)
        mask = self.table["transmission"] < min_thru
        # TODO: maybe use actually masked table here?
        self.table["transmission"][mask] = 0

    def fov_grid(self, which="waveset", **kwargs):
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release.", DeprecationWarning, stacklevel=2)
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = from_currsys(self.meta, self.cmds)
            # ..todo:: replace the 101 with a variable in !SIM
            wave = np.linspace(self.meta["wave_min"],
                               self.meta["wave_max"], 101)
            wave = quantify(wave, u.um)
            throughput = self.surface.transmission(wave)
            min_thru = self.meta["minimum_throughput"]
            valid_waves = np.where(throughput.value > min_thru)[0]
            if len(valid_waves) > 0:
                wave_edges = [min(wave[valid_waves].value),
                              max(wave[valid_waves].value)] * u.um
            else:
                raise ValueError("No transmission found above the threshold {}"
                                 " in this wavelength range {}. Did you open "
                                 "the shutter?"
                                 "".format(self.meta["minimum_throughput"],
                                           [self.meta["wave_min"],
                                            self.meta["wave_max"]]))
        else:
            wave_edges = []

        return wave_edges

    @property
    def fwhm(self):
        wave = self.surface.wavelength
        # noinspection PyProtectedMember
        thru = self.surface._get_ter_property("transmission", fmt="array")
        mask = thru >= 0.5
        if any(mask):
            dwave = wave[mask][-1] - wave[mask][0]
        else:
            dwave = 0 * wave.unit

        return dwave

    @property
    def centre(self):
        wave = self.surface.wavelength
        # noinspection PyProtectedMember
        thru = self.surface._get_ter_property("transmission", fmt="array")
        num = np.trapz(thru * wave**2, x=wave)
        den = np.trapz(thru * wave, x=wave)

        return num / den

    @property
    def center(self):
        return self.centre


class TopHatFilterCurve(FilterCurve):
    """
    A simple Top-Hat filter profile.

    Parameters
    ----------
    transmission : float
        [0..1] Peak transmission of filter

    blue_cutoff, red_cutoff : float
        [um] Blue and Red cutoff wavelengths

    wing_transmission : float, optional
        [0..1] Default 0. Wing transmission of filter outside the cutoff range

    Examples
    --------
    ::

        name: J_band_tophat
        class: TopHatFilterCurve
        kwargs:
            transmission : 0.9
            wing_transmission : 0.001
            blue_cutoff : 1.15
            red_cutoff : 1.35


    """

    required_keys = {"transmission", "blue_cutoff", "red_cutoff"}

    def __init__(self, cmds=None, **kwargs):
        check_keys(kwargs, self.required_keys, action="error")
        self.cmds = cmds

        wave_min = from_currsys("!SIM.spectral.wave_min", self.cmds)
        wave_max = from_currsys("!SIM.spectral.wave_max", self.cmds)
        blue = kwargs["blue_cutoff"]
        red = kwargs["red_cutoff"]
        peak = kwargs["transmission"]
        wing = kwargs.get("wing_transmission", 0)

        waveset = [wave_min, 0.999*blue, blue, red, red*1.001, wave_max]
        transmission = [wing, wing, peak, peak, wing, wing]

        tbl = Table(names=["wavelength", "transmission"],
                    data=[waveset, transmission])
        super().__init__(table=tbl, wavelength_unit="um",
                         action="transmission", cmds=self.cmds)
        self.meta.update(kwargs)


class DownloadableFilterCurve(FilterCurve):
    required_keys = {"filter_name", "filename_format"}

    def __init__(self, **kwargs):
        check_keys(kwargs, self.required_keys, action="error")
        filt_str = kwargs["filename_format"].format(kwargs["filter_name"])
        tbl = download_svo_filter(filt_str, return_style="table")
        super().__init__(table=tbl, **kwargs)


class SpanishVOFilterCurve(FilterCurve):
    """
    Pulls a filter transmission curve down from the Spanish VO filter service.

    Parameters
    ----------
    observatory : str
    instrument : str
    filter_name : str

    Examples
    --------
    ::

        name: HAWKI-Ks
        class: SpanishVOFilterCurve
        kwargs:
            observatory : Paranal
            instrument : HAWKI
            filter_name : Ks

    """

    required_keys = {"observatory", "instrument", "filter_name"}

    def __init__(self, **kwargs):
        check_keys(kwargs, self.required_keys, action="error")
        filt_str = "{}/{}.{}".format(kwargs["observatory"],
                                     kwargs["instrument"],
                                     kwargs["filter_name"])
        kwargs["name"] = kwargs["filter_name"]
        kwargs["svo_id"] = filt_str

        tbl = download_svo_filter(filt_str, return_style="table")
        super().__init__(table=tbl, **kwargs)


class FilterWheelBase(Effect):
    """Base class for Filter Wheels."""

    z_order: ClassVar[tuple[int, ...]] = (124, 224, 524)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = True
    report_table_rounding: ClassVar[int] = 4
    _current_str = "current_filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(kwargs, self.required_keys, action="error")

        self.meta.update(kwargs)

        self.filters = {}

    def apply_to(self, obj, **kwargs):
        """Use apply_to of current filter."""
        return self.current_filter.apply_to(obj, **kwargs)

    @property
    def surface(self):
        return self.current_filter.surface

    @property
    def throughput(self):
        return self.current_filter.throughput

    def fov_grid(self, which="waveset", **kwargs):
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release.", DeprecationWarning, stacklevel=2)
        return self.current_filter.fov_grid(which=which, **kwargs)

    def change_filter(self, filtername=None):
        """Change the current filter."""
        if filtername in self.filters.keys():
            self.meta["current_filter"] = filtername
        else:
            raise ValueError(f"Unknown filter requested: {filtername}")

    def add_filter(self, newfilter, name=None):
        """
        Add a filter to the FilterWheel.

        Parameters
        ----------
        newfilter : FilterCurve
        name : string
           Name to be used for the new filter. If `None` a name from
           the newfilter object is used.
        """
        if name is None:
            name = newfilter.display_name
        self.filters[name] = newfilter

    @property
    def current_filter(self):
        filter_eff = None
        filt_name = from_currsys(self.meta["current_filter"], self.cmds)
        if filt_name is not None:
            filter_eff = self.filters[filt_name]
        return filter_eff

    def __getattr__(self, item):
        return getattr(self.current_filter, item)

    def plot(self, which="x", wavelength=None, *, axes=None, **kwargs):
        """Plot TER curves.

        Parameters
        ----------
        which : {"x", "t", "e", "r"}, optional
            "x" plots throughput. "t","e","r" plot trans/emission/refl.
            Can be a combination, e.g. "tr" or "tex" to plot each.
        wavelength : array_like, optional
            DESCRIPTION. The default is None.
        axes : matplotlib axes, optional
            If given, plot into existing axes. The default is None.

        Returns
        -------
        fig : matplotlib figure
            Figure containing plots.

        """
        if axes is None:
            fig, axes = figure_factory(len(which), 1, iterable_axes=True)
        else:
            fig = axes.figure
            _guard_plot_axes(which, axes)

        for ter, ax in zip(which, axes):
            for name, _filter in self.filters.items():
                _filter.plot(which=ter, wavelength=wavelength, axes=ax,
                             plot_kwargs={"label": name}, **kwargs)

        fig.legend()
        return fig

    def get_table(self):
        names = list(self.filters.keys())
        ters = self.filters.values()
        centres = u.Quantity([ter.centre for ter in ters])
        widths = u.Quantity([ter.fwhm for ter in ters])
        blue = centres - 0.5 * widths
        red = centres + 0.5 * widths

        tbl = Table(names=["name", "centre", "width", "blue cutoff", "red cutoff"],
                    data=[names, centres, widths, blue, red])

        return tbl


class FilterWheel(FilterWheelBase):
    """
    Wheel holding a selection of predefined filters.

    Examples
    --------
    ::

        name: filter_wheel
        class: FilterWheel
        kwargs:
            filter_names: []
            filename_format: "filters/{}.
            current_filter: "Ks"

    """

    required_keys = {"filter_names", "filename_format", "current_filter"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        params = {"path": ""}
        self.meta.update(params)
        self.meta.update(kwargs)

        path = self._get_path()
        for name in from_currsys(self.meta["filter_names"], self.cmds):
            kwargs["name"] = name
            self.filters[name] = FilterCurve(filename=str(path).format(name),
                                             **kwargs)

        self.table = self.get_table()


class TopHatFilterWheel(FilterWheelBase):
    """
    A selection of top-hat filter curves as defined in the input lists.

    Parameters
    ----------
    filter_names: list of string

    transmissions: list of floats
        [0..1] Peak transmissions inside the cutoff limits

    wing_transmissions: list of floats
        [0..1] Wing transmissions outside the cutoff limits

    blue_cutoffs: list of floats
        [um]

    red_cutoffs: list of floats
        [um]

    current_filter: str, optional
        Name of current filter at initialisation. If no name is given, the
        first entry in `filter_names` is used by default.

    Examples
    --------
    ::

        name: top_hat_filter_wheel
        class: TopHatFilterWheel
        kwargs:
            filter_names: ["J", "H", "K"]
            transmissions: [0.9, 0.95, 0.85]
            wing_transmissions: [0., 0., 0.001]
            blue_cutoffs: [1.15, 1.45, 1.9]
            red_cutoffs: [1.35, 1.8, 2.4]
            current_filter: "K"

    """

    required_keys = {"filter_names", "transmissions", "wing_transmissions",
                     "blue_cutoffs", "red_cutoffs"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        current_filter = kwargs.get("current_filter",
                                    kwargs["filter_names"][0])
        params = {"current_filter": current_filter}
        self.meta.update(params)
        self.meta.update(kwargs)

        for i_filt, name in enumerate(self.meta["filter_names"]):
            effect_kwargs = {
                "name": name,
                "transmission": self.meta["transmissions"][i_filt],
                "wing_transmission": self.meta["wing_transmissions"][i_filt],
                "blue_cutoff": self.meta["blue_cutoffs"][i_filt],
                "red_cutoff": self.meta["red_cutoffs"][i_filt]}
            self.filters[name] = TopHatFilterCurve(**effect_kwargs)


class SpanishVOFilterWheel(FilterWheelBase):
    """
    A FilterWheel that loads all the filters from the Spanish VO service.

    .. warning::
       This use ``astropy.download_file(..., cache=True)``.

    The filter transmission curves probably won't change, but if you notice
    discrepancies, try clearing the astropy cache::

        >> from astropy.utils.data import clear_download_cache
        >> clear_download_cache()

    Parameters
    ----------
    observatory : str

    instrument : str

    current_filter : str
        Default filter name

    include_str, exclude_str : str
        String sequences that can be used to include or exclude filter names
        which contain a certain string.
        E.g. GTC/OSIRIS has curves for ``sdss_g`` and ``sdss_g_filter``.
        We can force the inclusion of only the filter curves by setting
        ``list_include_str: "_filter"``.

    Examples
    --------
    ::

        name: svo_filter_wheel
        class: SpanishVOFilterWheel
        kwargs:
            observatory: "GTC"
            instrument: "OSIRIS"
            current_filter: "sdss_r_filter"
            include_str: "_filter"

    """

    required_keys = {"observatory", "instrument", "current_filter"}

    def __init__(self, **kwargs):        
        super().__init__(**kwargs)

        params = {"include_str": None,         # passed to
                  "exclude_str": None,
                  }
        self.meta.update(params)
        self.meta.update(kwargs)

        obs, inst = self.meta["observatory"], self.meta["instrument"]
        filter_names = download_svo_filter_list(
            obs, inst, short_names=True,
            include=self.meta["include_str"], exclude=self.meta["exclude_str"])

        self.meta["filter_names"] = filter_names
        for name in filter_names:
            self.filters[name] = SpanishVOFilterCurve(observatory=obs,
                                                      instrument=inst,
                                                      filter_name=name)

        self.filters["open"] = FilterCurve(
            array_dict={"wavelength": [0.3, 3.0], "transmission": [1., 1.]},
            wavelength_unit="um", name="unity transmission")

        self.table = self.get_table()


class PupilTransmission(TERCurve):
    """
    Wavelength-independent transmission curve.

    Use this class to describe a cold stop or pupil mask that is
    characterised by "grey" transmissivity.
    The emissivity is set to zero, assuming that the mask is cold.
    """

    def __init__(self, transmission, cmds=None, **kwargs):
        self.params = {"wave_min": "!SIM.spectral.wave_min",
                       "wave_max": "!SIM.spectral.wave_max"}
        self.params.update(kwargs)
        self.cmds = cmds
        wave_min = from_currsys(self.params["wave_min"], self.cmds) * u.um
        wave_max = from_currsys(self.params["wave_max"], self.cmds) * u.um
        transmission = from_currsys(transmission, cmds=self.cmds)

        super().__init__(wavelength=[wave_min, wave_max],
                         transmission=[transmission, transmission],
                         emissivity=[0., 0.], **self.params)

    def update_transmission(self, transmission, **kwargs):
        self.__init__(transmission, **kwargs)


class ADCWheel(Effect):
    """
    Wheel holding a selection of predefined atmospheric dispersion correctors.

    Example
    -------
    ::

       name : adc_wheel
       class: ADCWheel
       kwargs:
           adc_names: []
           filename_format: "TER_ADC_{}.dat"
           current_adc: "const_90"
    """

    required_keys = {"adc_names", "filename_format", "current_adc"}
    z_order: ClassVar[tuple[int, ...]] = (125, 225, 525)
    report_plot_include: ClassVar[bool] = False
    report_table_include: ClassVar[bool] = True
    report_table_rounding: ClassVar[int] = 4
    _current_str = "current_adc"

    def __init__(self, cmds=None, **kwargs):
        super().__init__(cmds=cmds, **kwargs)
        check_keys(kwargs, self.required_keys, action="error")

        params = {"path": ""}
        self.meta.update(params)
        self.meta.update(kwargs)

        path = self._get_path()
        self.adcs = {}
        for name in from_currsys(self.meta["adc_names"], cmds=self.cmds):
            kwargs["name"] = name
            self.adcs[name] = TERCurve(filename=str(path).format(name),
                                       cmds=cmds,
                                       **kwargs)

        self.table = self.get_table()

    def apply_to(self, obj, **kwargs):
        """Use ``apply_to`` of current ADC."""
        return self.current_adc.apply_to(obj, **kwargs)

    def change_adc(self, adcname=None):
        """Change the current ADC."""
        if not adcname or adcname in self.adcs.keys():
            self.meta["current_adc"] = adcname
            self.include = adcname
        else:
            raise ValueError(f"Unknown ADC requested: {adcname}")

    @property
    def current_adc(self):
        """Return the currently used ADC."""
        curradc = from_currsys(self.meta["current_adc"], cmds=self.cmds)
        if not curradc:
            return False
        return self.adcs[curradc]

    def __getattr__(self, item):
        return getattr(self.current_adc, item)

    def get_table(self):
        """Create a table of ADCs with maximum throughput."""
        names = list(self.adcs.keys())
        adcs = self.adcs.values()
        tmax = np.array([adc.data["transmission"].max() for adc in adcs])

        tbl = Table(names=["name", "max_transmission"],
                    data=[names, tmax])
        return tbl


def _guard_plot_axes(which, axes):
    if len(which) > 1:
        if not isinstance(axes, Collection):
            raise TypeError(("axes must be collection of axes if which "
                             "contains more than one element"))
        if not len(axes) == len(which):
            raise ValueError("len of which and axes must match")
    else:
        if isinstance(axes, Collection):
            raise TypeError(("axes must be a single axes object if which "
                             "contains only one element"))
