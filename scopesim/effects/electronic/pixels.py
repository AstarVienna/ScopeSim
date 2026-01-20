# -*- coding: utf-8 -*-
"""Any effects regarding pixels, like binning and border effects."""

from typing import ClassVar

from .. import Effect
from ...detector import Detector
from ...utils import from_currsys, figure_factory, check_keys, real_colname
from .. import logger


class ReferencePixelBorder(Effect):
    """Remove signal from reference pixels.

    Detectors often have a number of rows and columns around the edges masked.
    These pixels serve as reference pixels for various purposes. They do not
    get any signal, but have all the detector effects, such as dark current
    and readout noise.

    .. versionchanged:: 0.11.2

       Re-implemented the effect with a new YAML syntax, see #840 for details.

    Parameters
    ----------
    border : list(4)
       a list with the number of rows and columns to be masked. The sequence
       should be [bottom, left, top, right]
    """

    z_order: ClassVar[tuple[int, ...]] = (861,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["border_sequence"] = "bottom left top right"
        if "border" not in self.meta:
            self.meta["border"] = [0, 0, 0, 0]
        else:
            self.meta["border"] = from_currsys(self.meta["border"], self.cmds)
        if isinstance(self.meta["border"], dict):
            for val in self.meta["border"].values():
                if len(val) != 4:
                    raise ValueError("All entries for 'border' must have exactly four values.")
        else:
            if len(self.meta["border"]) != 4:
                raise ValueError(
                    "Parameter 'border' must have exactly four values.")

    def apply_to(self, obj, **kwargs):
        """Mask border pixels."""
        if not isinstance(obj, Detector):
            logger.warning(
                "ReferencePixelBorder: got non-detector object: %s", type(obj))
            return obj

        logger.info(f"Applying border {from_currsys(self.meta['border'])}")
        if hasattr(self.meta["border"], "dic"):
            dtcr_id = obj.meta[real_colname("id", obj.meta)]
            border = self.meta["border"].dic[dtcr_id]
        elif isinstance(self.meta["border"], list):
            border = self.meta["border"]
        else:
            raise ValueError(
                "<ReferenceBorderPixel>.meta['border'] must be either "
                f"dict or list, but is {self.meta['border']}")

        if border[0] > 0:
            obj.data[:border[0], :] = 0
        if border[1] > 0:
            obj.data[:, :border[1]] = 0
        if border[2] > 0:
            obj.data[-border[2]:, :] = 0
        if border[3] > 0:
            obj.data[:, -border[3]:] = 0
        return obj

    def plot(self, det, **kwargs):
        """Show the masked detector image."""
        det = self.apply_to(det)
        _, ax = figure_factory()
        ax.imshow(det.data, origin="bottom", **kwargs)

    def __str__(self) -> str:
        """Return str(self)."""
        msg = (
            f"{self.__class__.__name__}: \"{self.display_name}\"\n"
            f"    {from_currsys(self.meta['border'], self.cmds)}"
            f"   ({self.meta['border_sequence']})\n"
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
