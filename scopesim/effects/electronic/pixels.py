# -*- coding: utf-8 -*-
"""Any effects regarding pixels, like binning and border effects."""

from typing import ClassVar

from .. import Effect
from ...base_classes import DetectorBase, ImagePlaneBase
from ...utils import from_currsys, figure_factory, check_keys


class ReferencePixelBorder(Effect):
    z_order: ClassVar[tuple[int, ...]] = (780,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        val = int(kwargs.get("all", 0))
        widths = {key: val for key in ["top", "bottom", "left", "right"]}
        self.meta.update(widths)
        self.meta.update(kwargs)

    def apply_to(self, implane, **kwargs):
        # .. todo: should this be ImagePlaneBase here?
        if isinstance(implane, ImagePlaneBase):
            if self.meta["top"] > 0:
                implane.hdu.data[:, -self.meta["top"]:] = 0
            if self.meta["bottom"] > 0:
                implane.hdu.data[:, :self.meta["bottom"]] = 0
            if self.meta["right"] > 0:
                implane.hdu.data[-self.meta["right"]:, :] = 0
            if self.meta["left"] > 0:
                implane.hdu.data[:self.meta["left"], :] = 0

        return implane

    def plot(self, implane, **kwargs):
        implane = self.apply_to(implane)
        fig, ax = figure_factory()
        ax.imshow(implane.data, origin="bottom", **kwargs)
        # fig.show()


class BinnedImage(Effect):
    required_keys = {"bin_size"}
    z_order: ClassVar[tuple[int, ...]] = (870,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
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
        if isinstance(det, DetectorBase):
            bx = from_currsys(self.meta["binx"], self.cmds)
            by = from_currsys(self.meta["biny"], self.cmds)
            image = det._hdu.data
            h, w = image.shape
            new_image = image.reshape((h//bx, bx, w//by, by))
            det._hdu.data = new_image.sum(axis=3).sum(axis=1)

        return det
