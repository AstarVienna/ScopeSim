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
            "action": "emission",
            "position": 0,  # position in surface table
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        self.compute_emission()


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
        tbl.add_column(np.zeros_like(lam).value, name="transmission")
        tbl.add_column(flux, name="emission")
        tbl.meta["wavelength_unit"] = tbl["wavelength"].unit
        tbl.meta["emission_unit"] = tbl["emission"].unit

        self.surface.table = tbl
        self.surface.meta.update(tbl.meta)
