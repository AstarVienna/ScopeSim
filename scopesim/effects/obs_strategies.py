# -*- coding: utf-8 -*-
"""
Effects describing observing strategies.

- ChopNodCombiner: simulate chop-nod cycle
"""

from typing import ClassVar

import numpy as np

from scopesim.base_classes import DetectorBase
from scopesim.effects import Effect
from scopesim.utils import from_currsys, check_keys


class ChopNodCombiner(Effect):
    """
    Creates and combines 4 images for each of the chop/nod positions.

    - AA : original position ``(dx, dy) = (0, 0)``
    - AB : chop position ``(dx, dy) = chop_offsets``
    - BA : nod position ``(dx, dy) = nod_offsets``
    - BB : chop-nod position ``(dx, dy) = nod_offsets + chop_offsets``

    Images are combined using::

        im_combined = (AA - AB) - (BA - BB)

    If no ``nod_offset`` is given, it is set to the inverse of ``chop_offset``.

    ``ChopNodCombiner`` is a detector effect and should be placed last in the
     detector yaml (after the noise effects).

    Keyword Arguments
    -----------------
    chop_offsets : tuple, optional
        [arcsec] (dx, dy) offset of chop position relative to AA
    nod_offsets : tuple, optional
        [arcsec] (dx, dy) offset of nod position relative to AA

    Example yaml entry
    ------------------
    ::

        name: perpendicular_chop_nod_slanted_pattern
        description: chop throw to (+5, 0) and nod throw to (-7, +10) arcsec
        class: ChopNodCombiner
        include: True
        kwargs:
            pixel_scale : "!INST.pixel_scale"
            chop_offsets : (5, 0)
            nod_offsets : (-7, 10)

    """

    required_keys = {"chop_offsets", "pixel_scale"}
    z_order: ClassVar[tuple[int, ...]] = (863,)

    def __init__(self, **kwargs):
        check_keys(kwargs, self.required_keys)

        super().__init__(**kwargs)
        params = {
            "chop_offsets": None,
            "nod_offsets": None,
            "pixel_scale": None,
            "include": True,
        }
        self.meta.update(params)
        self.meta.update(kwargs)

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, DetectorBase):
            return obj

        chop_offsets = from_currsys(self.meta["chop_offsets"], self.cmds)
        nod_offsets = from_currsys(self.meta["nod_offsets"], self.cmds)
        if nod_offsets is None:
            nod_offsets = -np.array(chop_offsets)

        # these offsets are in pixels, not in arcsec or mm
        pixel_scale = float(from_currsys(self.meta["pixel_scale"], self.cmds))
        chop_offsets_px = (np.array(chop_offsets) / pixel_scale).astype(int)
        nod_offsets_px = (np.array(nod_offsets) / pixel_scale).astype(int)

        image = obj.hdu.data
        obj.hdu.data = chop_nod_image(image, chop_offsets_px, nod_offsets_px)

        return obj


def chop_nod_image(img, chop_offsets, nod_offsets=None):
    """Create four copies and combine in chop-nod pattern."""
    if nod_offsets is None:
        nod_offsets = tuple(-np.array(chop_offsets))

    im_aa = np.copy(img)
    im_ab = np.roll(im_aa, chop_offsets, (1, 0))
    im_ba = np.roll(im_aa, nod_offsets, (1, 0))
    im_bb = np.roll(im_ba, chop_offsets, (1, 0))

    im_comb = (im_aa - im_ab) - (im_ba - im_bb)
    return im_comb
