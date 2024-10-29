# -*- coding: utf-8 -*-
"""
Effects related to rotation of the field/CCD.

Classes:
- RotateCCD - Rotates CCD by integer multiples of 90 degrees
"""

from typing import ClassVar

import numpy as np

from . import Effect
from ..utils import from_currsys, check_keys
from ..base_classes import DetectorBase


class Rotate90CCD(Effect):
    """
    Rotates CCD by integer multiples of 90 degrees.

    rotations kwarg is number of counter-clockwise rotations

    Author: Dave jones

    """

    required_keys = {"rotations"}
    z_order: ClassVar[tuple[int, ...]] = (809,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        """See parent docstring."""
        if isinstance(obj, DetectorBase):
            rotations = from_currsys(self.meta["rotations"], self.cmds)
            obj._hdu.data = np.rot90(obj._hdu.data, rotations)

        return obj
