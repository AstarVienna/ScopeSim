"""Transmission, emissivity, reflection curves"""
import logging
from pathlib import Path

import numpy as np
import skycalc_ipy
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

from .effects import Effect
from .ter_curves_utils import add_edge_zeros
from .ter_curves_utils import combine_two_spectra, apply_throughput_to_cube
from .ter_curves_utils import download_svo_filter, download_svo_filter_list
from ..base_classes import SourceBase, FOVSetupBase
from ..optics.surface import SpectralSurface
from ..source.source import Source
from ..utils import from_currsys, quantify, check_keys, find_file


class TERCurve(Effect):
    """
    Transmission, Emissivity, Reflection Curve

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
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)
        params = {"z_order": [10, 110, 510],
                  "ignore_wings": False,
                  "wave_min": "!SIM.spectral.wave_min",
                  "wave_max": "!SIM.spectral.wave_max",
                  "wave_unit": "!SIM.spectral.wave_unit",
                  "wave_bin": "!SIM.spectral.spectral_bin_width",
                  "bg_cell_width": "!SIM.computing.bg_cell_width",
                  "report_plot_include": True,
                  "report_table_include": False}
        self.meta.update(params)
        self.meta.update(kwargs)

        self.surface = SpectralSurface()
        self.surface.meta.update(self.meta)
        self._background_source = None

        data = self.get_data()
        if self.meta["ignore_wings"]:
            data = add_edge_zeros(data, "wavelength")
        if data is not None:
            # Assert that get_data() did not give us an image.
            assert isinstance(data, Table), "TER Curves must be tables."
            self.surface.table = data
            self.surface.table.meta.update(self.meta)

    # ####### added in new branch

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, SourceBase):
            assert isinstance(obj, Source), "Only Source supported."
            self.meta = from_currsys(self.meta)
            wave_min = quantify(self.meta["wave_min"], u.um).to(u.AA)
            wave_max = quantify(self.meta["wave_max"], u.um).to(u.AA)

            thru = self.throughput

            # apply transmission to source spectra
            for isp, spec in enumerate(obj.spectra):
                obj.spectra[isp] = combine_two_spectra(spec, thru, "multiply",
                                                       wave_min, wave_max)

            # apply transmission to cube fields
            for icube, cube in enumerate(obj.cube_fields):
                obj.cube_fields[icube] = apply_throughput_to_cube(cube, thru)

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
            # add a single pixel ImageHDU for the extended background with a
            # size of 1 degree
            # bg_cell_width = from_currsys(self.meta["bg_cell_width"])

            flux = self.emission
            bg_hdu = fits.ImageHDU()
            # TODO: The make_imagehdu_from_table below has been replaced with
            #       the empty ImageHDU above in fbca416. That change might,
            #       have been fine (or not?), but now there is no use anywhere
            #       in the code of make_imagehdu_from_table or bg_cell_width,
            #       so maybe these need to be removed?
            # bg_hdu = make_imagehdu_from_table([0], [0], [1], bg_cell_width * u.arcsec)

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

    # #######

    def plot(self, which="x", wavelength=None, ax=None, new_figure=True,
             label=None, **kwargs):
        """

        Parameters
        ----------
        which : str
            "x" plots throughput. "t","e","r" plot trans/emission/refl
        wavelength : list, np.ndarray
        ax : matplotlib.Axis
        new_figure : start a new figure (or add to the existing one)
        label : the label to use (ignored)
        kwargs

        Returns
        -------

        """
        import matplotlib.pyplot as plt
        if new_figure:
            # from matplotlib import rcParams
            # rcParams["axes.formatter.useoffset"] = False
            plt.figure(figsize=(10, 5))

        self.meta.update(kwargs)
        params = from_currsys(self.meta)

        for ii, ter in enumerate(which):
            if ax is None:
                plt.subplot(len(which), 1, ii+1)

            if wavelength is None:
                wunit = params["wave_unit"]
                wave = np.arange(quantify(params["wave_min"], wunit).value,
                                 quantify(params["wave_max"], wunit).value,
                                 quantify(params["wave_bin"], wunit).value)
                wave *= u.Unit(wunit)
            else:
                wave = wavelength

            plot_kwargs = self.meta.get("plot_kwargs", {})
            surf = self.surface

            if "t" in ter:
                y = surf.transmission(wave)
            elif "e" in ter:
                y = surf.emission(wave)
            elif "r" in ter:
                y = surf.reflection(wave)
            else:
                y = surf.throughput(wave)

            plt.plot(wave, y, **plot_kwargs)

            wave_unit = self.meta.get("wavelength_unit")
            plt.xlabel(f"Wavelength [{wave_unit}]")
            y_str = {"t": "Transmission", "e": "Emission",
                     "r": "Reflectivity", "x": "Throughput"}
            plt.ylabel(f"{y_str[ter]} [{y.unit}]")

        return plt.gcf()


class AtmosphericTERCurve(TERCurve):
    def __init__(self, **kwargs):
        super(AtmosphericTERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [111, 511]
        self.meta["action"] = "transmission"
        self.meta["position"] = 0       # position in surface table
        self.meta.update(kwargs)
        self.surface.meta.update(self.meta)


class SkycalcTERCurve(AtmosphericTERCurve):
    def __init__(self, **kwargs):
        """
        Retrieves an atmospheric spectrum from ESO's skycalc server

        kwarg parameters
        ----------------
        skycalc parameters can be found by calling::

            >>> import skycalc_ipy
            >>> skycalc_ipy.SkyCalc().keys

        .. note:: Different to skycalc_ipy, wmin and wmax must be given in units
            of ``um``

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
        super(SkycalcTERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [112, 512]
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
        return from_currsys(self.meta["include"])

    @include.setter
    def include(self, item):
        self.meta["include"] = item
        if item is True and self.skycalc_table is None:
            self.load_skycalc_table()

    def load_skycalc_table(self):
        use_local_file = from_currsys(self.meta["use_local_skycalc_file"])
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
            scale_factor = u.Unit(from_currsys(self.meta["wunit"])).to(u.nm)
            for key in ["wmin", "wmax", "wdelta"]:
                if key in self.meta:
                    self.meta[key] = from_currsys(self.meta[key]) * scale_factor

        conn_kwargs = {key: self.meta[key] for key in self.meta
                       if key in self.skycalc_conn.defaults}
        conn_kwargs = from_currsys(conn_kwargs)
        self.skycalc_conn.values.update(conn_kwargs)

        try:
            tbl = self.skycalc_conn.get_sky_spectrum(return_type="table")
        except ConnectionError:
            msg = "Could not connect to skycalc server"
            logging.exception(msg)
            raise ValueError(msg)

        return tbl


class QuantumEfficiencyCurve(TERCurve):
    def __init__(self, **kwargs):
        super(QuantumEfficiencyCurve, self).__init__(**kwargs)
        self.meta["action"] = "transmission"
        self.meta["z_order"] = [113, 513]
        self.meta["position"] = -1          # position in surface table


class FilterCurve(TERCurve):
    """
    Parameters
    ----------
    position : int, optional
    filter_name : str, optional
        ``Ks`` - corresponding to the filter name in the filename pattern
    filename_format : str, optional
        ``TC_filter_{}.dat``

    Can either be created using the standard 3 options:
    - ``filename``: direct filename of the filter curve
    - ``table``: an ``astropy.Table``
    - ``array_dict``: a dictionary version of a table: ``{col_name1: values, }``

    or by passing the combination of ``filter_name`` and ``filename_format`` as
    kwargs. Here all filter file names follow a pattern (e.g. see above) and the
    ``{}`` are replaced by ``filter_name`` at run time. ``filter_name`` can
    also be a !bang string for a ``__currsys__`` entry: ``"!INST.filter_name"``

    """
    def __init__(self, **kwargs):
        if not np.any([key in kwargs for key in ["filename", "table",
                                                 "array_dict"]]):
            if "filter_name" in kwargs and "filename_format" in kwargs:
                filt_name = from_currsys(kwargs["filter_name"])
                file_format = from_currsys(kwargs["filename_format"])
                kwargs["filename"] = file_format.format(filt_name)
            else:
                raise ValueError("FilterCurve must be passed one of (`filename`"
                                 " `array_dict`, `table`) or both "
                                 f"(`filter_name`, `filename_format`): {kwargs}")

        super(FilterCurve, self).__init__(**kwargs)
        if self.table is None:
            raise ValueError("Could not initialise filter. Either filename not "
                             "found, or array are not compatible")

        params = {"minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "action": "transmission",
                  "position": -1,               # position in surface table
                  "wing_flux_level": None,
                  "name": "untitled filter"}
        self.meta.update(params)
        self.meta["z_order"] = [114, 214, 514]
        self.meta.update(kwargs)

        min_thru = from_currsys(self.meta["minimum_throughput"])
        mask = self.table["transmission"] < min_thru
        self.table["transmission"][mask] = 0

    def fov_grid(self, which="waveset", **kwargs):
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = from_currsys(self.meta)
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
                raise ValueError("No transmission found above the threshold {} "
                                 "in this wavelength range {}. Did you open "
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
    A simple Top-Hat filter profile

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
    def __init__(self, **kwargs):
        required_keys = ["transmission", "blue_cutoff", "red_cutoff"]
        check_keys(kwargs, required_keys, action="error")

        wave_min = from_currsys("!SIM.spectral.wave_min")
        wave_max = from_currsys("!SIM.spectral.wave_max")
        blue = kwargs["blue_cutoff"]
        red = kwargs["red_cutoff"]
        peak = kwargs["transmission"]
        wing = kwargs.get("wing_transmission", 0)

        waveset = [wave_min, 0.999*blue, blue, red, red*1.001, wave_max]
        transmission = [wing, wing, peak, peak, wing, wing]

        tbl = Table(names=["wavelength", "transmission"],
                    data=[waveset, transmission])
        super(TopHatFilterCurve, self).__init__(table=tbl,
                                                wavelength_unit="um",
                                                action="transmission")
        self.meta.update(kwargs)


class DownloadableFilterCurve(FilterCurve):
    def __init__(self, **kwargs):
        required_keys = ["filter_name", "filename_format"]
        check_keys(kwargs, required_keys, action="error")
        filt_str = kwargs["filename_format"].format(kwargs["filter_name"])
        tbl = download_svo_filter(filt_str, return_style="table")
        super(DownloadableFilterCurve, self).__init__(table=tbl, **kwargs)


class SpanishVOFilterCurve(FilterCurve):
    """
    Pulls a filter transmission curve down from the Spanish VO filter service

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
    def __init__(self, **kwargs):
        required_keys = ["observatory", "instrument", "filter_name"]
        check_keys(kwargs, required_keys, action="error")
        filt_str = "{}/{}.{}".format(kwargs["observatory"],
                                     kwargs["instrument"],
                                     kwargs["filter_name"])
        kwargs["name"] = kwargs["filter_name"]
        kwargs["svo_id"] = filt_str

        tbl = download_svo_filter(filt_str, return_style="table")
        super(SpanishVOFilterCurve, self).__init__(table=tbl, **kwargs)


class FilterWheel(Effect):
    """
    This wheel holds a selection of predefined filters.

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

    def __init__(self, **kwargs):
        required_keys = ["filter_names", "filename_format", "current_filter"]
        check_keys(kwargs, required_keys, action="error")

        super(FilterWheel, self).__init__(**kwargs)

        params = {"z_order": [124, 224, 524],
                  "path": "",
                  "report_plot_include": True,
                  "report_table_include": True,
                  "report_table_rounding": 4}
        self.meta.update(params)
        self.meta.update(kwargs)

        path = Path(self.meta["path"], from_currsys(self.meta["filename_format"]))
        self.filters = {}
        for name in from_currsys(self.meta["filter_names"]):
            kwargs["name"] = name
            self.filters[name] = FilterCurve(filename=str(path).format(name),
                                             **kwargs)

        self.table = self.get_table()

    def apply_to(self, obj, **kwargs):
        """Use apply_to of current filter"""
        return self.current_filter.apply_to(obj, **kwargs)

    def fov_grid(self, which="waveset", **kwargs):
        return self.current_filter.fov_grid(which=which, **kwargs)

    def change_filter(self, filtername=None):
        """Change the current filter"""
        if filtername in self.filters.keys():
            self.meta["current_filter"] = filtername
        else:
            raise ValueError(f"Unknown filter requested: {filtername}")

    def add_filter(self, newfilter, name=None):
        """
        Add a filter to the FilterWheel

        Parameters
        ==========
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
        filt_name = from_currsys(self.meta["current_filter"])
        if filt_name is not None:
            filter_eff = self.filters[filt_name]
        return filter_eff

    @property
    def display_name(self):
        return (f"{self.meta['name']} : "
                f"[{from_currsys(self.meta['current_filter'])}]")

    def __getattr__(self, item):
        return getattr(self.current_filter, item)

    def plot(self, which="x", wavelength=None, **kwargs):
        """

        Parameters
        ----------
        which : str
            "x" plots throughput. "t","e","r" plot trans/emission/refl
        wavelength
        kwargs

        Returns
        -------

        """
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 5))

        for ii, ter in enumerate(which):
            ax = plt.subplot(len(which), 1, ii+1)
            for name, _filter in self.filters.items():
                _filter.plot(which=ter, wavelength=wavelength, ax=ax,
                             new_figure=False, plot_kwargs={"label": name},
                             **kwargs)

        # plt.semilogy()
        plt.legend()

        return plt.gcf()

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


class TopHatFilterWheel(FilterWheel):
    """
    A selection of top-hat filter curves as defined in the input lists

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
        first entry in ``filter_names`` is used by default.

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
    def __init__(self, **kwargs):
        required_keys = ["filter_names", "transmissions", "wing_transmissions",
                         "blue_cutoffs", "red_cutoffs"]
        check_keys(kwargs, required_keys, action="error")

        super(FilterWheel, self).__init__(**kwargs)

        params = {"z_order": [124, 224, 524],
                  "report_plot_include": True,
                  "report_table_include": True,
                  "report_table_rounding": 4,
                  "current_filter": kwargs.get("current_filter",
                                               kwargs["filter_names"][0])
                  }
        self.meta.update(params)
        self.meta.update(kwargs)

        self.filters = {}
        for i in range(len(self.meta["filter_names"])):
            name = self.meta["filter_names"][i]
            effect_kwargs = {"name": name,
                             "transmission": self.meta["transmissions"][i],
                             "wing_transmission": self.meta["wing_transmissions"][i],
                             "blue_cutoff": self.meta["blue_cutoffs"][i],
                             "red_cutoff": self.meta["red_cutoffs"][i]}
            self.filters[name] = TopHatFilterCurve(**effect_kwargs)


class SpanishVOFilterWheel(FilterWheel):
    """
    A FilterWheel that loads all the filters from the Spanish VO service

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
    def __init__(self, **kwargs):
        required_keys = ["observatory", "instrument", "current_filter"]
        check_keys(kwargs, required_keys, action="error")

        # Call Effect.init, *NOT* FilterWheel.init --> different required_keys
        super(FilterWheel, self).__init__(**kwargs)

        params = {"z_order": [124, 224, 524],
                  "report_plot_include": True,
                  "report_table_include": True,
                  "report_table_rounding": 4,
                  "include_str": None,         # passed to
                  "exclude_str": None,
                  }
        self.meta.update(params)
        self.meta.update(kwargs)

        obs, inst = self.meta["observatory"], self.meta["instrument"]
        inc, exc = self.meta["include_str"], self.meta["exclude_str"]
        filter_names = download_svo_filter_list(obs, inst, short_names=True,
                                                include=inc, exclude=exc)

        self.meta["filter_names"] = filter_names
        self.filters = {name: SpanishVOFilterCurve(observatory=obs,
                                                   instrument=inst,
                                                   filter_name=name)
                        for name in filter_names}
        self.filters["open"] = FilterCurve(array_dict={"wavelength": [0.3, 3.0],
                                                       "transmission": [1., 1.]},
                                           wavelength_unit="um",
                                           name="unity transmission")

        self.table = self.get_table()


class PupilTransmission(TERCurve):
    """
    Wavelength-independent transmission curve

    Use this class to describe a cold stop or pupil mask that is
    characterised by "grey" transmissivity.
    The emissivity is set to zero, assuming that the mask is cold.
    """
    def __init__(self, transmission, **kwargs):
        self.params = {"wave_min": "!SIM.spectral.wave_min",
                       "wave_max": "!SIM.spectral.wave_max"}
        self.params.update(kwargs)
        wave_min = from_currsys(self.params["wave_min"]) * u.um
        wave_max = from_currsys(self.params["wave_max"]) * u.um
        transmission = from_currsys(transmission)

        super().__init__(wavelength=[wave_min, wave_max],
                         transmission=[transmission, transmission],
                         emissivity=[0., 0.], **self.params)

    def update_transmission(self, transmission, **kwargs):
        self.__init__(transmission, **kwargs)


class ADCWheel(Effect):
    """
    This wheel holds a selection of predefined atmospheric dispersion
    correctors.

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
    def __init__(self, **kwargs):
        required_keys = ["adc_names", "filename_format", "current_adc"]
        check_keys(kwargs, required_keys, action="error")

        super().__init__(**kwargs)

        params = {"z_order": [125, 225, 525],
                  "path": "",
                  "report_plot_include": False,
                  "report_table_include": True,
                  "report_table_rounding": 4}
        self.meta.update(params)
        self.meta.update(kwargs)

        path = Path(self.meta["path"], from_currsys(self.meta["filename_format"]))
        self.adcs = {}
        for name in from_currsys(self.meta["adc_names"]):
            kwargs["name"] = name
            self.adcs[name] = TERCurve(filename=str(path).format(name),
                                       **kwargs)

        self.table = self.get_table()

    def apply_to(self, obj, **kwargs):
        """Use apply_to of current adc"""
        return self.current_adc.apply_to(obj, **kwargs)

    def change_adc(self, adcname=None):
        """Change the current ADC"""
        if not adcname or adcname in self.adcs.keys():
            self.meta["current_adc"] = adcname
            self.include = adcname
        else:
            raise ValueError(f"Unknown ADC requested: {adcname}")

    @property
    def current_adc(self):
        """Return the currently used ADC"""
        curradc = from_currsys(self.meta["current_adc"])
        if not curradc:
            return False
        return self.adcs[curradc]

    @property
    def display_name(self):
        return (f"{self.meta['name']} : "
                f"[{from_currsys(self.meta['current_adc'])}]")

    def __getattr__(self, item):
        return getattr(self.current_adc, item)

    def get_table(self):
        """Create a table of ADCs with maximum throughput"""
        names = list(self.adcs.keys())
        adcs = self.adcs.values()
        tmax = np.array([adc.data["transmission"].max() for adc in adcs])

        tbl = Table(names=["name", "max_transmission"],
                    data=[names, tmax])
        return tbl
