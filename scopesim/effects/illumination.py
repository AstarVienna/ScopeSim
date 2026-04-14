# -*- coding: utf-8 -*-
"""Image-plane illumination effects"""

from typing import Callable, ClassVar

import numpy as np
from astropy.modeling.functional_models import Gaussian2D

from . import Effect
from ..optics.image_plane import ImagePlane
from ..utils import figure_factory


__all__ = ["Illumination", "gaussian2d", "quadratic_vignetting"]


def gaussian2d(shape, amp=1.0, mu=(0.0, 0.0), sigma=(2000.0, 2000.0), theta=0.0):
    """Normalised 2D elliptical Gaussian.

    Parameters
    ----------
    shape : tuple of int
        (ny, nx) image shape in pixels.
    amp : float
        Peak amplitude (normalised to 1 by default).
    mu : tuple of float
        (x, y) centre offset in pixels from image centre.
    sigma : tuple of float
        (sx, sy) Gaussian widths in pixels.
    theta : float
        Rotation angle in radians, counterclockwise.
    """
    ny, nx = shape[-2], shape[-1]
    y, x = np.ogrid[:ny, :nx]
    x = x - nx / 2
    y = y - ny / 2
    model = Gaussian2D(amplitude=amp, x_mean=mu[0], y_mean=mu[1],
                       x_stddev=sigma[0], y_stddev=sigma[1], theta=theta)
    return model(x, y)


def quadratic_vignetting(shape, falloff=0.01, r_ref=None, mu=(0.0, 0.0),
                    stretch=(1.0, 1.0, 1.0, 1.0)):
    """Quadratic vignetting pattern with independent stretch factors.

    Parameters
    ----------
    shape : tuple of int
        (ny, nx) image shape in pixels.
    max_falloff : float
        Fractional illumination drop at ``r_ref`` (e.g. 0.01 = 1 %).
    r_ref : float or None
        Reference radius in stretched pixels.  Defaults to the corner distance.
    mu : tuple of float
        (x, y) offset of the vignetting centre in pixels from the image centre.
    stretch : tuple of float
        ``(+x, -x, +y, -y)`` independent scale factors for
        half-planes respectively.  All 1.0 gives a
        circular pattern.  A value > 1 widens the falloff in that direction
        (shallower); < 1 narrows it (steeper).
    """
    ny, nx = shape[-2], shape[-1]
    
    yy, xx = np.ogrid[:ny, :nx]
    dx = xx - (nx / 2 + mu[0])
    dy = yy - (ny / 2 + mu[1])

    sx = np.where(dx >= 0, stretch[0], stretch[1])
    sy = np.where(dy >= 0, stretch[2], stretch[3])

    r2 = (dx / sx)**2 + (dy / sy)**2

    if r_ref is None:
        r2_ref = r2.max()
    else:
        r2_ref = r_ref**2

    return np.clip(1.0 - falloff * r2 / r2_ref, 0.0, 1.0)


class Illumination(Effect):
    """Large-scale illumination variation across the image plane.

    Parameters
    ----------
    model : callable, optional
        Function ``f(shape, **kwargs) -> ndarray`` returning the
        illumination map. Defaults to :func:`gaussian2d`.
    modelargs : dict, optional
        Keyword arguments forwarded to ``model``. If omitted, the model's
        own defaults are used.

    include : str
        Turn effect on/off from the IRDB
        default.yaml.  Defaults to ``"!DET.include_illumination"``.

    Examples
    --------

    Polynomial vignetting with <1 % falloff (auto r_ref from image shape)::

        eff = Illumination(model=poly_vignetting, modelargs={"max_falloff": 0.01})

    Custom model::

        def my_model(shape, slope=-0.001):
            ny, nx = shape[-2], shape[-1]
            y, x = np.ogrid[:ny, :nx]
            r = np.sqrt((x - nx / 2)**2 + (y - ny / 2)**2)
            return np.clip(1 + slope * r, 0, None)

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
        illumination_map = self._model(shape, **self._modelargs)
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
