# -*- coding: utf-8 -*-
"""Image-plane illumination effects"""

from typing import Callable, ClassVar

import numpy as np

from . import Effect
from ..optics.image_plane import ImagePlane
from ..utils import figure_factory


__all__ = ["Illumination", "gaussian2d", "poly_vignetting"]


def gaussian2d(xx, yy, amp=1.0, mu=(0.0, 0.0), sigma=(2000.0, 2000.0)):
    """Normalised 2D elliptical Gaussian.

    Parameters
    ----------
    xx, yy : ndarray
        2D coordinate grids in pixels relative to the image centre.
    amp : float
        Peak amplitude (normalised to 1 by default).
    mu : tuple of float
        (x, y) centre offset in pixels from image centre.
    sigma : tuple of float
        (sx, sy) Gaussian widths in pixels.
    """
    dx = xx - mu[0]
    dy = yy - mu[1]
    return amp * np.exp(-(dx**2 / (2 * sigma[0]**2) + dy**2 / (2 * sigma[1]**2)))


def poly_vignetting(xx, yy, max_falloff=0.01, r_ref=None, mu=(0.0, 0.0),
                    stretch=(1.0, 1.0, 1.0, 1.0)):
    """Asymmetric polynomial vignetting with direct falloff control.

    Parameters
    ----------
    xx, yy : ndarray
        2D coordinate grids in pixels relative to the image centre.
    max_falloff : float
        Fractional illumination drop at ``r_ref`` (e.g. 0.01 = 1 %).
    r_ref : float or None
        Reference radius in stretched pixels.  Defaults to the corner distance.
    mu : tuple of float
        (x, y) offset of the vignetting centre in pixels from the image centre.
    stretch : tuple of float
        ``(sx_pos, sx_neg, sy_pos, sy_neg)`` — independent scale factors for
        the +x, -x, +y, and -y half-planes respectively.  All 1.0 gives a
        circular pattern.  A value > 1 widens the falloff in that direction
        (shallower); < 1 narrows it (steeper).
    """
    dx = xx - mu[0]
    dy = yy - mu[1]
    sx = np.where(dx >= 0, stretch[0], stretch[1])
    sy = np.where(dy >= 0, stretch[2], stretch[3])
    r2 = (dx / sx)**2 + (dy / sy)**2
    if r_ref is None:
        r_ref = np.sqrt(r2.max())
    return np.clip(1.0 - max_falloff * r2 / r_ref**2, 0.0, 1.0)


class Illumination(Effect):
    """Large-scale illumination variation across the image plane.

    Parameters
    ----------
    model : callable, optional
        Function ``f(xx, yy, **kwargs) -> ndarray`` returning the
        illumination map. Defaults to :func:`gaussian2d`.
    modelargs : dict, optional
        Keyword arguments forwarded to ``model``. If omitted, the model's
        own defaults are used.

    include : str
        Bang-string reference to toggle the effect on/off from the IRDB
        default.yaml.  Defaults to ``"!DET.include_illumination"``.

    Examples
    --------
    IRDB default.yaml entry::

        - name: illumination
          description: Large-scale illumination variation
          class: Illumination
          kwargs:
            model: gaussian2d
            modelargs:
              sigma: [2000, 2000]
            include: "!DET.include_illumination"

    Polynomial vignetting with <1 % falloff (auto r_ref from image shape)::

        eff = Illumination(model=poly_vignetting, modelargs={"max_falloff": 0.01})

    Custom model::

        def my_model(xx, yy, slope=-0.001):
            return np.clip(1 + slope * np.sqrt(xx**2 + yy**2), 0, None)

        eff = Illumination(model=my_model, modelargs={"slope": -0.0005})
    """

    z_order: ClassVar[tuple[int, ...]] = (750,)

    def __init__(self, model: Callable = gaussian2d, modelargs: dict = None, **kwargs):
        super().__init__(**kwargs)
        self.meta.setdefault("include", "!DET.include_illumination")
        self._model = model
        self._modelargs = modelargs or {}
        self._map = None
        self._map_shape = None

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, ImagePlane):
            return obj

        shape = obj.hdu.data.shape
        if self._map is None or shape != self._map_shape:
            self._map = self._make_map(shape)
            self._map_shape = shape

        obj.hdu.data *= self._map
        return obj

    def _make_map(self, shape):
        ny, nx = shape[-2], shape[-1]
        xx, yy = np.meshgrid(np.arange(nx) - nx / 2, np.arange(ny) - ny / 2)
        illumination_map = self._model(xx, yy, **self._modelargs)
        return illumination_map.astype(np.float32)

    def plot(self):
        if self._map is None:
            raise RuntimeError("No illumination map cached — run a simulation first.")

        fig, ax = figure_factory()
        im = ax.imshow(self._map, origin="lower", vmin=0.98, vmax=1, cmap="gray_r")
        fig.colorbar(im, ax=ax, label="Relative illumination")
        ax.set_title("Illumination")
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")
        return fig
