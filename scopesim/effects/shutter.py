# -*- coding: utf-8 -*-
"""Contains the Shutter effect."""

from typing import ClassVar

from . import Effect
from ..base_classes import ImagePlaneBase
from ..utils import get_logger


logger = get_logger(__name__)


class Shutter(Effect):
    """Simulate a closed shutter, useful for dark exposures."""

    z_order: ClassVar[tuple[int, ...]] = (799,)

    def apply_to(self, obj, **kwargs):
        """Set all pixels of image plane to zero."""
        if not isinstance(obj, ImagePlaneBase):
            return obj

        logger.warning("Shutter is closed, setting all pixels to zero.")
        obj.data[:] = 0.0
        return obj
