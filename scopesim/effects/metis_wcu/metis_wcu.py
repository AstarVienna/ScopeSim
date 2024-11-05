"""Classes for the METIS Warm Calibration Unit"""

import numpy as np
from ..ter_curves import TERCurve
from ...utils import get_logger

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
            "z_order": [113, 513]
        }
        self.meta.update(params)
        self.meta.update(kwargs)
