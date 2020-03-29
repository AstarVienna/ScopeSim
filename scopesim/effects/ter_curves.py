import numpy as np
from astropy import units as u
from astropy.table import Table
from os import path as pth

from synphot import SourceSpectrum

from .effects import Effect
from ..optics.surface import SpectralSurface
from ..utils import from_currsys, quantify, check_keys
from .ter_curves_utils import download_svo_filter


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
        self.meta["z_order"] = [10, 110]
        self.meta["ignore_wings"] = False
        self.meta.update(kwargs)
        self.surface = SpectralSurface()
        self.surface.meta.update(self.meta)
        data = self.get_data()
        if self.meta["ignore_wings"]:
            data = add_edge_zeros(data, "wavelength")
        if data is not None:
            self.surface.table = data
            self.surface.table.meta.update(self.meta)


def add_edge_zeros(tbl, wave_colname):
    if isinstance(tbl, Table):
        vals = np.zeros(len(tbl.colnames))
        col_i = np.where(col == wave_colname for col in tbl.colnames)[0][0]
        sgn = np.sign(np.diff(tbl[wave_colname][:2]))
        vals[col_i] = tbl[wave_colname][0] * (1 - 1e-7 * sgn)
        tbl.insert_row(0, vals)
        vals[col_i] = tbl[wave_colname][-1] * (1 + 1e-7 * sgn)
        tbl.insert_row(len(tbl), vals)

    return tbl


class AtmosphericTERCurve(TERCurve):
    def __init__(self, **kwargs):
        super(AtmosphericTERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [111]
        self.meta["action"] = "transmission"
        self.meta["area"] = "!TEL.area"
        self.meta["area_unit"] = "m2"
        self.meta["position"] = 0       # position in surface table
        self.meta.update(kwargs)


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
        self.meta["z_order"] = [112]
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
        tbl = self.skycalc_conn.get_sky_spectrum(return_type="table",
                                                 filename=filename)
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
        self.meta["z_order"] = [113]
        self.meta["position"] = -1          # position in surface table


class FilterCurve(TERCurve):
    """
    kwargs
    ------
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
        params = {"minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "action": "transmission",
                  "position": -1,               # position in surface table
                  "wing_flux_level": None}
        self.meta.update(params)
        self.meta["z_order"] = [114, 214]
        self.meta.update(kwargs)

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

