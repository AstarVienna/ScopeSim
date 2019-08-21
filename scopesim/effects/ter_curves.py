import numpy as np
from astropy import units as u

import skycalc_ipy

from .effects import Effect
from ..optics.surface import SpectralSurface
from ..utils import from_currsys, quantify


class TERCurve(Effect):
    """
    Transmission, Emissivity, Reflection Curve

    Must contain a wavelength column, and one or more of the following:
    ``transmission``, ``emissivity``, ``reflection``.
    Additionally in the header there
    should be the following keywords: wavelength_unit

    """
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [10, 210]
        self.surface = SpectralSurface()
        self.surface.meta.update(self.meta)
        data = self.get_data()
        if data is not None:
            self.surface.table = data


class AtmosphericTERCurve(TERCurve):
    def __init__(self, **kwargs):
        super(AtmosphericTERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [11, 211]
        self.meta["area"] = "!TEL.area"
        self.meta["area_unit"] = "m2"
        self.meta.update(kwargs)


class SkycalcTERCurve(TERCurve):
    def __init__(self, **kwargs):
        """
        Retrieves an atmospheric spectrum from ESO's skycalc server

        kwarg parameters
        ----------------
        skycalc parameters can be found by calling::

            >>> skycalc_ipy.SkyCalc().keys

        .. note:: Compared to skycalc_ipy, wmin and wmax must be given in units
            of ``um``

        """

        super(SkycalcTERCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [12, 212]
        self.meta["action"] = "transmission"
        self.meta["area"] = "!TEL.area"
        self.meta["area_unit"] = "m2"
        self.meta.update(kwargs)

        self.skycalc_conn = skycalc_ipy.SkyCalc()
        self.query_server(**kwargs)

        if "name" not in self.meta:
            self.meta["name"] = self.skycalc_conn["observatory"]

    def query_server(self, **kwargs):
        kwargs = from_currsys(kwargs)
        self.skycalc_conn.values.update(kwargs)

        tbl = self.skycalc_conn.get_sky_spectrum(return_type="table")
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
        self.meta["z_order"] = [13, 213]


class FilterCurve(TERCurve):
    def __init__(self, **kwargs):
        if np.all([key not in kwargs for key in ["filename", "table",
                                                 "array_dict"]]):
            if "filter_name" in kwargs and "filename_format" in kwargs:
                filt_name = from_currsys(kwargs["filter_name"])
                file_format = kwargs["filename_format"]
                kwargs["filename"] = file_format.format(filt_name)
            else:
                raise ValueError("FilterCurve must be passed one of (`filename`"
                                 " `array_dict`, `table`) or both "
                                 "(`filter_name`, `filename_format`):"
                                 "{}".format(kwargs))

        super(FilterCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [14, 214]
        self.meta["min_throughput"] = "!SIM.spectral.minimum_throughput"
        self.meta["action"] = "transmission"

    def fov_grid(self, which="waveset", **kwargs):
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = from_currsys(self.meta)
            # ..todo:: replace the 101 with a variable in !SIM
            wave = np.linspace(self.meta["wave_min"],
                               self.meta["wave_max"], 101)
            wave = quantify(wave, u.um).to(u.um)
            throughput = self.surface.transmission(wave)
            min_thru = self.meta["min_throughput"]
            valid_waves = np.where(throughput.value > min_thru)[0]
            if len(valid_waves) > 0:
                wave_edges = [min(wave[valid_waves].value),
                              max(wave[valid_waves].value)] * u.um
            else:
                raise ValueError("No transmission found above the threshold {} "
                                 "in this wavelength range {}. Did you open "
                                 "the shutter?"
                                 "".format(self.meta["min_throughput"],
                                           [self.meta["wave_min"],
                                            self.meta["wave_max"]]))
        else:
            wave_edges = []

        return wave_edges
