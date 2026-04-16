# -*- coding: utf-8 -*-
"""Image-plane illumination effects."""

from typing import ClassVar
from collections.abc import Callable, Mapping

import numpy as np
from astropy import units as u
from astropy.modeling.functional_models import Gaussian2D

from . import Effect
from ..optics.image_plane import ImagePlane
from ..utils import figure_factory


__all__ = ["Illumination", "gaussian2d", "quadratic_vignetting"]


def gaussian2d(
    shape: tuple[int, int],
    amp: float = 1.0,
    mu: tuple[float, float] = (0.0, 0.0),
    sigma: tuple[float, float] = (2000.0, 2000.0),
    theta: u.Quantity[u.deg] | float = 0.0 * u.deg,
) -> np.ndarray:
    """
    Normalised 2D elliptical Gaussian to be used for vignetting map.

    .. versionadded:: PLACEHOLDER_NEXT_RELEASE_VERSION

    Parameters
    ----------
    shape : tuple[int, int]
        Image shape in pixels (ny, nx).
    amp : float, optional
        Peak amplitude. The default is 1.0 (normalized).
    mu : tuple[float, float], optional
        Offset of the peak center in pixels (x, y) from the image center.
        The default is (0.0, 0.0), i.e. no offset.
    sigma : tuple[float, float], optional
        Gaussian widths in pixels (sx, sy). The default is (2000.0, 2000.0).
    theta : float, optional
        Rotation angle (if float, the angle is interpreted in degrees),
        counterclockwise. The default is 0°.

    Returns
    -------
    np.ndarray
        Vignetting map.

    """
    nx, ny = reversed(shape)
    y, x = np.ogrid[:ny, :nx]
    x = x - nx / 2
    y = y - ny / 2

    model = Gaussian2D(
        amplitude=amp,
        x_mean=mu[0],
        y_mean=mu[1],
        x_stddev=sigma[0],
        y_stddev=sigma[1],
        theta=theta << u.deg,
    )
    return model(x, y)


def quadratic_vignetting(
    shape: tuple[int, int],
    falloff: float = 0.01,
    r_ref: float | None = None,
    mu: tuple[float, float] = (0.0, 0.0),
    stretch: tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0),
) -> np.ndarray:
    """
    Quadratic vignetting pattern with independent stretch factors.

    .. versionadded:: PLACEHOLDER_NEXT_RELEASE_VERSION

    Parameters
    ----------
    shape : tuple[int, int]
        Image shape in pixels (ny, nx).
    falloff : float, optional
        Fractional illumination drop at `r_ref`. The default is 0.01 (= 1 %).
    r_ref : float | None, optional
        Reference radius in stretched pixels. If None (the default), use the
        corner distance.
    mu : tuple[float, float], optional
        Offset of the vignetting center in pixels (x, y) from the image center.
        The default is (0.0, 0.0), i.e. no offset.
    stretch : tuple[float, float, float, float], optional
        ``(+x, -x, +y, -y)`` independent scale factors for half-planes
        respectively. All 1.0 gives a circular pattern. A value > 1 widens the
        falloff in that direction (shallower); < 1 narrows it (steeper).
        The default is (1.0, 1.0, 1.0, 1.0).

    Returns
    -------
    np.ndarray
        Vignetting map.

    """
    nx, ny = reversed(shape)

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

    .. versionadded:: PLACEHOLDER_NEXT_RELEASE_VERSION

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
    Polynomial vignetting with <1 % falloff (auto r_ref from image shape)

    >>> eff = Illumination(model=poly_vignetting, modelargs={"falloff": 0.01})

    Custom model

    >>> def my_model(shape, slope=-0.001):
    >>>     ny, nx = shape[-2], shape[-1]
    >>>     y, x = np.ogrid[:ny, :nx]
    >>>     r = np.sqrt((x - nx / 2)**2 + (y - ny / 2)**2)
    >>>     return np.clip(1 + slope * r, 0, None)
    >>>
    >>> eff = Illumination(model=my_model, modelargs={"slope": -0.0005})

    """

    z_order: ClassVar[tuple[int, ...]] = (750,)

    def __init__(
        self,
        model: Callable = gaussian2d,
        modelargs: Mapping | None = None,
        **kwargs,
    ) -> None:
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
        """Plot effect."""
        if self._map is None:
            raise RuntimeError(
                "No illumination map cached — run a simulation first."
            )

        fig, ax = figure_factory()
        im = ax.imshow(
            self._map, origin="lower", vmin=0.98, vmax=1., cmap="gray_r",
        )
        fig.colorbar(im, ax=ax, label="Relative illumination")
        ax.set_title("Illumination")
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")
        return fig
