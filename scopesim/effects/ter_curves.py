import skycalc_ipy

from .effects import Effect
from ..optics.surface import SpectralSurface
from ..utils import from_currsys


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
