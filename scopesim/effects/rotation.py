"""Effects related to rotation of the field/CCD

Classes:
- RotateCCD - Rotates CCD by integer multiples of 90 degrees
"""

import logging
import numpy as np

from . import Effect
from .. import utils
from ..utils import from_currsys
from ..base_classes import DetectorBase


class Rotate90CCD(Effect):
    """
    Rotates CCD by integer multiples of 90 degrees
    rotations kwarg is number of counter-clockwise rotations
    
    Author: Dave jones

    """
    def __init__(self, **kwargs):
        super(Rotate90CCD, self).__init__(**kwargs)
        params = {"z_order": [809]}
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ["rotations"]
        utils.check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            rotations = from_currsys(self.meta["rotations"])
            obj._hdu.data = np.rot90(obj._hdu.data,rotations)

        return obj
