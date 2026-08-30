# -*- coding: utf-8 -*-
"""Any effects regarding pixels, like binning and border effects."""

from typing import ClassVar

from numpy.typing import ArrayLike, NDArray

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

    def __call__(self, data: ArrayLike, border: dict[int, int]) -> NDArray:
        if border[0] > 0:
            data[:border[0], :] = 0
        if border[1] > 0:
            data[:, :border[1]] = 0
        if border[2] > 0:
            data[-border[2]:, :] = 0
        if border[3] > 0:
            data[:, -border[3]:] = 0
        return data

    def _get_border(self, det_meta) -> dict[int, int]:
        if hasattr(self.meta["border"], "dic"):
            dtcr_id = det_meta[real_colname("id", det_meta)]
            return self.meta["border"].dic[dtcr_id]
        if isinstance(self.meta["border"], list):
            return self.meta["border"]
        raise ValueError(
            f"{self.__class__.__name__}.meta['border'] must be either "
            f"dict or list, but is {self.meta['border']}")

    def apply_to(self, obj, **kwargs):
        """Mask border pixels."""
        if not isinstance(obj, Detector):
            logger.warning(
                "ReferencePixelBorder: got non-detector object: %s", type(obj))
            return obj

        logger.info(f"Applying border {from_currsys(self.meta['border'])}")
        obj.data = self(obj.data, self._get_border(obj.meta))
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


class BinnedImageBase(Effect):
    """Base class for binning effects."""

    z_order: ClassVar[tuple[int, ...]] = (870,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def __call__(self, data: ArrayLike) -> NDArray:
        raise NotImplementedError("Subclasses must implement this.")

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            return det

        det._hdu.data = self(det._hdu.data)
        return det


class BinnedImage(BinnedImageBase):
    required_keys = {"bin_size"}

    def __call__(self, data: ArrayLike) -> NDArray:
        height, width = data.shape
        binned_data = data.reshape((
            height//self.bin_size,
            self.bin_size,
            width//self.bin_size,
            self.bin_size
        )).sum(axis=3).sum(axis=1)
        return binned_data

    @property
    def bin_size(self) -> int:
        return from_currsys(self.meta["bin_size"], self.cmds)


class UnequalBinnedImage(BinnedImageBase):
    required_keys = {"binx","biny"}

    def __call__(self, data: ArrayLike) -> NDArray:
        height, width = data.shape
        binned_data = data.reshape((
            height//self.binx,
            self.binx,
            width//self.biny,
            self.binx
        )).sum(axis=3).sum(axis=1)
        return binned_data

    @property
    def binx(self) -> int:
        return from_currsys(self.meta["binx"], self.cmds)

    @property
    def biny(self) -> int:
        return from_currsys(self.meta["biny"], self.cmds)
