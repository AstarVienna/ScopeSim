import numpy as np
from astropy import units as u
from os import path as pth
import warnings

from astropy.io import fits
from astropy.table import Table
from astropy import units as u

from synphot import SourceSpectrum
from synphot.units import PHOTLAM

from .ter_curves_utils import combine_two_spectra, add_edge_zeros, download_svo_filter
from .effects import Effect
from ..optics.surface import SpectralSurface
from ..source.source_utils import make_imagehdu_from_table
from ..source.source import Source
from ..base_classes import SourceBase
from ..utils import from_currsys, quantify, check_keys


class TERCurve(Effect):
    """
    Transmission, Emissivity, Reflection Curve

    Must contain a wavelength column, and one or more of the following:
    ``transmission``, ``emissivity``, ``reflection``.
    Additionally in the header there
    should be the following keywords: wavelength_unit

    kwargs that can be passed::

        "rescale_emission" : { "filter_name": str, "value": float, "unit": str}

    Examples
    --------
    Inside a YAML file description::

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
                filter_name: "Paranal/HAWKI.Ks"
                value: 15.5
                unit: ABmag

    Inside an ASCII file::

        # name: bogus_surface
        # wavelength_unit: um
        # emission_unit: ph s-1 m-2 um-1
        # rescale_emission: {filter_name: "Paranal/HAWKI.Ks", value: 36.3, unit: Jy}
        wavelength  transmission    emission
        0.3         0.9             1
        3.0         0.9             1

    """
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)
        params = {"z_order": [10, 110, 510],
                  "ignore_wings": False,
                  "wave_min": "!SIM.spectral.wave_min",
                  "wave_max": "!SIM.spectral.wave_max",
                  "wave_unit": "!SIM.spectral.wave_unit",
                  "wave_bin": "!SIM.spectral.spectral_resolution",
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
            self.surface.table = data
            self.surface.table.meta.update(self.meta)

    # ####### added in new branch

    def apply_to(self, obj):
        if isinstance(obj, SourceBase):
            self.meta = from_currsys(self.meta)
            wave_min = quantify(self.meta["wave_min"], u.um).to(u.AA)
            wave_max = quantify(self.meta["wave_max"], u.um).to(u.AA)

            # apply transmission to source spectra
            for ii, spec in enumerate(obj.spectra):
                thru = self.throughput
                obj.spectra[ii] = combine_two_spectra(spec, thru, "multiply",
                                                      wave_min, wave_max)

            # add the effect background to the source background field
            if self.background_source is not None:
                obj.append(self.background_source)

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
            flux = self.emission
            bg_hdu = make_imagehdu_from_table([0], [0], [1], 1 * u.deg)
            bg_hdu.header["BG_SRC"] = True
            bg_hdu.header["BG_SURF"] = self.meta.get("name", "<untitled>")
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
            plt.xlabel("Wavelength [{}]".format(wave_unit))
            y_str = {"t": "Transmission", "e": "Emission",
                     "r": "Reflectivity", "x": "Throughput"}
            plt.ylabel("{} [{}]".format(y_str[ter], y.unit))

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

        .. note:: Compared to skycalc_ipy, wmin and wmax must be given in units
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
        import skycalc_ipy

        super(SkycalcTERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [112, 512]
        self.meta.update(kwargs)

        self.skycalc_conn = skycalc_ipy.SkyCalc()
        self.query_server()
        if "name" not in self.meta:
            self.meta["name"] = self.skycalc_conn["observatory"]

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

        local_path = from_currsys("!SIM.file.local_packages_path")
        filename = pth.join(local_path, "skycalc_temp.fits")
        try:
            tbl = self.skycalc_conn.get_sky_spectrum(return_type="table",
                                                     filename=filename)
        except:
            warnings.warn("Could not connect to skycalc server")
            if pth.exists(filename):
                pass
            else:
                raise ValueError("No local copy exists: {}".format(filename))

        for i, colname in enumerate(["wavelength", "transmission", "emission"]):
            tbl.columns[i].name = colname
        tbl.meta["wavelength_unit"] = tbl.columns[0].unit
        tbl.meta["emission_unit"] = tbl.columns[2].unit
        self.surface.table = tbl
        self.surface.meta.update(tbl.meta)


class QuantumEfficiencyCurve(TERCurve):
    def __init__(self, **kwargs):
        super(QuantumEfficiencyCurve, self).__init__(**kwargs)
        self.meta["action"] = "transmission"
        self.meta["z_order"] = [113, 513]
        self.meta["position"] = -1          # position in surface table


class FilterCurve(TERCurve):
    """
    Other Parameters
    ----------------
    position : int
    filter_name : str
        ``Ks`` - corresponding to the filter name in the filename pattern
    filename_format : str
        ``TC_filter_{}.dat``

    Can either be created using the standard 3 options:
    - ``filename``: direct filename of the filer curve
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
                                 "(`filter_name`, `filename_format`):"
                                 "{}".format(kwargs))

        super(FilterCurve, self).__init__(**kwargs)
        if self.table is None:
            raise ValueError("Could not initialise filter. Either filename not "
                             "found, or array are not compatible")

        params = {"minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "action": "transmission",
                  "position": -1,               # position in surface table
                  "wing_flux_level": None}
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
        thru = self.surface._get_ter_property("transmission", fmt="array")
        num = np.trapz(thru * wave**2, x=wave)
        den = np.trapz(thru * wave, x=wave)

        return num / den

    @property
    def center(self):
        return self.centre


class DownloadableFilterCurve(FilterCurve):
    def __init__(self, **kwargs):
        required_keys = ["filter_name", "filename_format"]
        check_keys(kwargs, required_keys, action="error")
        filt_str = kwargs["filename_format"].format(kwargs["filter_name"])
        tbl = download_svo_filter(filt_str, return_style="table")
        super(DownloadableFilterCurve, self).__init__(table=tbl, **kwargs)


class SpanishVOFilterCurve(FilterCurve):
    def __init__(self, **kwargs):
        required_keys = ["observatory", "instrument", "filter_name"]
        check_keys(kwargs, required_keys, action="error")
        filt_str = "{}/{}.{}".format(kwargs["observatory"],
                                     kwargs["instrument"],
                                     kwargs["filter_name"])
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

        path = pth.join(self.meta["path"],
                        from_currsys(self.meta["filename_format"]))
        self.filters = {}
        for name in self.meta["filter_names"]:
            kwargs["name"] = name
            self.filters[name] = FilterCurve(filename=path.format(name), **kwargs)

        self.table = self.get_table()


    def apply_to(self, obj):
        '''Use apply_to of current filter'''
        return self.current_filter.apply_to(obj)

    def fov_grid(self, which="waveset", **kwargs):
        return self.current_filter.fov_grid(which=which, **kwargs)

    def change_filter(self, filtername=None):
        '''Change the current filter'''
        if filtername in self.filters.keys():
            self.meta['current_filter'] = filtername
        else:
            raise ValueError("Unknown filter requested: " + filtername)

    @property
    def current_filter(self):
        return self.filters[from_currsys(self.meta["current_filter"])]

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
            for name in self.filters:
                self.filters[name].plot(which=ter, wavelength=wavelength,
                                        ax=ax, new_figure=False,
                                        plot_kwargs={"label": name}, **kwargs)

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
        self.make_ter_curve(transmission)

    def update_transmission(self, transmission, **kwargs):
        self.params.update(kwargs)
        self.make_ter_curve(transmission)

    def make_ter_curve(self, transmission):
        wave_min = from_currsys(self.params["wave_min"]) * u.um
        wave_max = from_currsys(self.params["wave_max"]) * u.um
        transmission = from_currsys(transmission)

        super().__init__(wavelength=[wave_min, wave_max],
                         transmission=[transmission, transmission],
                         emissivity=[0., 0.], **self.params)

