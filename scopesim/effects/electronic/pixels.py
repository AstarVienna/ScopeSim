# -*- coding: utf-8 -*-
"""Any effects regarding pixels, like binning and border effects."""

from typing import ClassVar

from .. import Effect
from ...detector import Detector
from ...utils import from_currsys, figure_factory, check_keys
from .. import logger

class ReferencePixelBorder(Effect):
    """Remove signal from reference pixels

    Detectors often have a number of rows and columns around the edges masked.
    These pixels serve as reference pixels for various purposes. They do not
    get any signal, but have all the detector effects, such as dark current
    and readout noise.

    Parameters
    ----------
    border : list(4)
       a list with the number of rows and columns to be masked. The sequence
       should be [bottom, left, top, right]
    """
    z_order: ClassVar[tuple[int, ...]] = (780,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["bottom"] = 0
        self.meta["left"] = 0
        self.meta["top"] = 0
        self.meta["right"] = 0
        self.meta.update(kwargs)

        if "border" in self.meta:
            if len(self.meta["border"]) != 4:
                raise ValueError("Parameter 'border' must have exactly four entries.")
            self.meta['bottom'] = int(self.meta['border'][0])
            self.meta['left'] = int(self.meta['border'][1])
            self.meta['top'] = int(self.meta['border'][2])
            self.meta['right'] = int(self.meta['border'][3])

    def apply_to(self, obj, **kwargs):
        """Mask border pixels"""
        if not isinstance(obj, Detector):
            logger.warning("ReferencePixelBorder: got non-detector object")
            return obj

        if self.meta['bottom'] > 0:
            obj.data[:self.meta['bottom'], :] = 0
        if self.meta['left'] > 0:
            obj.data[:, :self.meta['left']] = 0
        if self.meta['top'] > 0:
            obj.data[-self.meta['top']:, :] = 0
        if self.meta['right'] > 0:
            obj.data[:, -self.meta['right']:] = 0
        return obj

    def plot(self, det, **kwargs):
        """Show the masked detector image"""
        det = self.apply_to(det)
        _, ax = figure_factory()
        ax.imshow(det.data, origin="bottom", **kwargs)

    def __str__(self) -> str:
        """Return str(self)."""
        msg = (
            f"{self.__class__.__name__}: \"{self.display_name}\"\n"
            f"   bottom:    {self.meta['bottom']}\n"
            f"   left:      {self.meta['left']}\n"
            f"   top:       {self.meta['top']}\n"
            f"   right:     {self.meta['right']}\n"
        )
        return msg



class BinnedImage(Effect):
    required_keys = {"bin_size"}
    z_order: ClassVar[tuple[int, ...]] = (870,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            return det

        bs = from_currsys(self.meta["bin_size"], self.cmds)
        image = det._hdu.data
        h, w = image.shape
        new_image = image.reshape((h//bs, bs, w//bs, bs))
        det._hdu.data = new_image.sum(axis=3).sum(axis=1)

        return det


class UnequalBinnedImage(Effect):
    required_keys = {"binx","biny"}
    z_order: ClassVar[tuple[int, ...]] = (870,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            return det

        bx = from_currsys(self.meta["binx"], self.cmds)
        by = from_currsys(self.meta["biny"], self.cmds)
        image = det._hdu.data
        h, w = image.shape
        new_image = image.reshape((h//bx, bx, w//by, by))
        det._hdu.data = new_image.sum(axis=3).sum(axis=1)

        return det
