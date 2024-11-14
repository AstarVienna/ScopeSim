"""Classes for the METIS Warm Calibration Unit"""

import numpy as np
from astropy.table import Table
from astropy.modeling.models import BlackBody
from astropy import units as u
from ..ter_curves import TERCurve
from ...utils import get_logger, seq

logger = get_logger(__name__)

class BlackBodySource(TERCurve):
    """Black Body Source

    This class returns a TERCurve that describes the output emission
    from the integrating sphere of the METIS WCU. The calculations include
    - black-body emission from the black-body source
    - coupling into the integrating sphere
    - integrating sphere magnification factor
    - thermal emission from the integrating sphere
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "z_order": [113, 513],
            "action": "emissivity",
            "position": 0,  # position in surface table
        }
        self.meta.update(params)
        self.meta.update(kwargs)
        self._background_source = None

        self.compute_emission()

    @property
    def emission(self):
        return self.surface.emission

    def set_temperature(self, bb_temp: [float | u.Quantity]=None,
                        wcu_temp: [float | u.Quantity]=None):
        """Change the black-body temperature

        Parameters
        ----------
        bb_temp, wcu_temp : float, Quantity
            new temperatures for the BB source and the ambient WCU, respectively.
            If float, the unit is assumed to be Kelvin.
        """
        if isinstance(bb_temp, float):
            bb_temp = bb_temp << u.K
            self.meta["bb_temp"] = bb_temp
        elif isinstance(bb_temp, u.Quantity):
            try:
                bb_temp = bb_temp.to(u.K, equivalencies=u.temperature())
                self.meta["bb_temp"] = bb_temp
            except u.UnitConversionError:
                logger.warning(
                    "UnitConversionError: Parameter bb_temp cannot be converted to Kelvin")
                return None

        if isinstance(wcu_temp, float):
            wcu_temp = wcu_temp << u.K
            self.meta["wcu_temp"] = wcu_temp
        elif isinstance(wcu_temp, u.Quantity):
            try:
                wcu_temp = wcu_temp.to(u.K, equivalencies=u.temperature())
                self.meta["wcu_temp"] = wcu_temp
            except u.UnitConversionError:
                logger.warning(
                    "UnitConversionError: Parameter wcu_temp cannot be converted to Kelvin")
                return None

        self.compute_emission()

        return None

    def compute_emission(self):
        """Compute the emission at the exit of the integrating sphere"""
        bb_temp = self.meta["bb_temp"] << u.K
        wcu_temp = self.meta["wcu_temp"] << u.K
        tbl = Table()
        lam = seq(2.2, 15, 0.01) * u.um
        bb_scale = 1 * u.ph / (u.s * u.m**2 * u.sr * u.um)
        bb_lam = (BlackBody(bb_temp, scale=bb_scale) +
                  BlackBody(wcu_temp, scale=bb_scale))
        flux = bb_lam(lam)
        tbl.add_column(lam, name="wavelength")
        tbl.add_column(np.ones_like(lam).value, name="transmission")
        tbl.add_column(flux, name="emission")
        tbl.meta["wavelength_unit"] = tbl["wavelength"].unit
        tbl.meta["emission_unit"] = tbl["emission"].unit

        self.surface.table = tbl
        self.surface.meta.update(tbl.meta)

    def __str__(self) -> str:
        return f"""{self.__class__.__name__}: "{self.display_name}"
        BlackBody temperature: {self.meta['bb_temp']}
        WCU temperature:       {self.meta['wcu_temp']}"""
