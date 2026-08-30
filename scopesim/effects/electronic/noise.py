# -*- coding: utf-8 -*-
"""Any kinds of electronic or photonic noise."""

from typing import ClassVar
from collections.abc import Mapping
from numbers import Real  # matches int, float and all the numpy scalars

import numpy as np
from numpy.typing import ArrayLike, NDArray
from astropy.io import fits

from .. import Effect
from ...detector import Detector
from ...utils import from_currsys, figure_factory, check_keys, real_colname
from . import logger


class Bias(Effect):
    """Adds a constant bias level to readout."""

    required_keys = {"bias"}
    z_order: ClassVar[tuple[int, ...]] = (855,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def __call__(self, data: ArrayLike) -> NDArray:
        return data + self.bias_level

    @property
    def bias_level(self) -> float:
        return from_currsys(self.meta["bias"], self.cmds)

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, Detector):
            return obj

        obj.data = self(obj.data)
        return obj


class PoorMansHxRGReadoutNoise(Effect):
    required_keys = {"noise_std", "n_channels", "ndit"}
    z_order: ClassVar[tuple[int, ...]] = (811,)
    report_plot_include: ClassVar[bool] = False
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "pedestal_fraction": 0.3,
            "read_fraction": 0.4,
            "line_fraction": 0.25,
            "channel_fraction": 0.05,
            "random_seed": "!SIM.random.seed",
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def __call__(self, data: ArrayLike, ndit: int = 1) -> NDArray:
        self.meta = from_currsys(self.meta, self.cmds)

        ron_frame = self._make_ron_frame(data.shape)

        stacked_ron_frame = np.zeros_like(ron_frame)
        for i in range(ndit):
            dx = np.random.randint(0, ron_frame.shape[1])
            dy = np.random.randint(0, ron_frame.shape[0])
            stacked_ron_frame += np.roll(ron_frame, (dy, dx), axis=(0, 1))

        # TODO: this .T is ugly. Work out where things are getting switched and remove it!
        return data + stacked_ron_frame.T

    @property
    def noise_std(self) -> int:
        return from_currsys(self.meta["noise_std"], self.cmds)

    @property
    def n_channels(self) -> int:
        return from_currsys(self.meta["n_channels"], self.cmds)

    @property
    def random_seed(self) -> int:
        return from_currsys(self.meta["random_seed"], self.cmds)

    @staticmethod
    def _pseudo_random_field(
        rng,
        scale: float = 1.,
        size: tuple[int, ...] = (1024, 1024),
    ) -> NDArray:
        n = 256
        image = np.zeros(size)
        batch = rng.normal(loc=0, scale=scale, size=(2*n, 2*n))
        for y in range(0, size[1], n):
            for x in range(0, size[0], n):
                i, j = rng.integers(n, size=2)
                dx, dy = min(size[0]-x, n), min(size[1]-y, n)
                image[x:x+dx, y:y+dy] = batch[i:i+dx, j:j+dy]

        return image

    def _make_ron_frame(self, shape: tuple[int, ...]) -> NDArray:
        # TODO: Add yaml-settable seed here
        rng = np.random.default_rng(self.random_seed)

        channel_fraction = self.meta["channel_fraction"]
        line_fraction = self.meta["line_fraction"]
        pedestal_fraction = self.meta["pedestal_fraction"]
        read_fraction = self.meta["read_fraction"]

        pixel_std = self.noise_std * (pedestal_fraction + read_fraction)**0.5
        if shape < (1024, 1024):
            pixel = rng.normal(loc=0, scale=pixel_std, size=shape)
            line = rng.normal(
                loc=0,
                scale=self.noise_std * line_fraction**0.5,
                size=shape[1],
            )
        else:
            # TODO: Why bother with this pseudo random function?
            pixel = self._pseudo_random_field(rng, scale=pixel_std, size=shape)
            line = pixel[0]

        channel = np.repeat(
            rng.normal(
                loc=0,
                scale=self.noise_std * channel_fraction**0.5,
                size=self.n_channels,
            ),
            max(1, shape[0] // self.n_channels) + 1,
            axis=0,
        )

        return (pixel + line).T + channel[:shape[0]]

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            return det

        ndit = from_currsys(self.meta["ndit"], self.cmds)
        det._hdu.data = self(det._hdu.data, ndit)
        return det

    def plot(self, det, **kwargs):
        """Plot effect image."""
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data, origin="lower")

    def plot_hist(self, det, **kwargs):
        """Plot effect histogram."""
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


class BasicReadoutNoise(Effect):
    """Readout noise computed as: ron * sqrt(NDIT)."""

    required_keys = {"noise_std", "ndit"}
    z_order: ClassVar[tuple[int, ...]] = (811,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def __call__(self, data: ArrayLike, ndit: int = 1) -> NDArray:
        rng = np.random.default_rng(self.random_seed)
        scale = self.noise_std * np.sqrt(float(ndit))
        return data + rng.normal(loc=0, scale=scale, size=data.shape)

    @property
    def noise_std(self) -> float:
        return from_currsys(self.meta["noise_std"], self.cmds)

    @property
    def random_seed(self) -> int:
        return from_currsys(self.meta["random_seed"], self.cmds)

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            return det

        ndit = from_currsys(self.meta["ndit"], self.cmds)
        det._hdu.data = self(det._hdu.data, ndit)
        return det

    def plot(self, det):
        """Plot effect image."""
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        """Plot effect histogram."""
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


# TODO: Is this really a "noise" effect? Sounds more like "electrons" tbh.
class PixelResponseNonUniformity(Effect):
    """Pixel Response Non-Uniformity (PRNU).

    Models the fixed pattern of per-pixel gain variations across the detector
    arising from manufacturing differences in quantum efficiency. Each pixel is
    multiplied by a gain factor drawn from N(1, ``prnu_std``) keyed by detector
    ID. The gain map is generated once per detector on first use and reused
    identically across all subsequent exposures.

    .. versionadded:: 0.11.3

    Parameters
    ----------
    prnu_std : float or dict
        Standard deviation of the per-pixel gain distribution.

    prnu_seed : int, fixed

    include:  "!DET.include_prnu"

    Example
    -------
    ::

       - name: prnu
         description: Pixel response non-uniformity
         class: PixelResponseNonUniformity
         kwargs:
           prnu_std: 0.001
           prnu_seed: 42
           include: "!DET.include_prnu"

    """

    required_keys: ClassVar[set] = set()
    z_order: ClassVar[tuple[int, ...]] = (805,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)
        self._gain_maps = {}  # keyed by dtcr_id

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, Detector):
            return obj

        random_seed = from_currsys(self.meta.get("prnu_seed"), self.cmds)
        id_key = real_colname("id", obj.meta)
        dtcr_id = obj.meta[id_key] if id_key is not None else None

        prnu_std_meta = from_currsys(self.meta["prnu_std"], self.cmds)
        if isinstance(prnu_std_meta, Mapping):
            prnu_std = float(from_currsys(prnu_std_meta[dtcr_id], self.cmds))
        elif isinstance(prnu_std_meta, Real):
            prnu_std = float(prnu_std_meta)
        else:
            raise TypeError(
                "<PixelResponseNonUniformity>.meta['prnu_std'] must be a float "
                f"or a dict keyed by detector ID, got {type(prnu_std_meta)}"
            )

        shape = obj.data.shape
        if dtcr_id not in self._gain_maps:
            rng = np.random.default_rng(random_seed)
            self._gain_maps[dtcr_id] = rng.normal(
                loc=1.0, scale=prnu_std, size=shape,
            )

        if self._gain_maps[dtcr_id].shape != shape:
            raise ValueError("gain map shape mismatch")

        obj.data = obj.data * self._gain_maps[dtcr_id]
        return obj

    def plot(self, det_id=None):
        """Plot effect."""
        if not self._gain_maps:
            raise RuntimeError("No gain map yet - run a simulation first.")
        key = det_id if det_id in self._gain_maps else next(iter(self._gain_maps))
        gain_map = self._gain_maps[key]
        dev = np.max(np.abs(gain_map - 1.0))
        fig, ax = figure_factory()
        im = ax.imshow(gain_map, origin="lower", aspect="auto",
                       vmin=1 - dev, vmax=1 + dev)
        fig.colorbar(im, ax=ax, label="per-pixel gain")
        return fig


class ShotNoise(Effect):
    z_order: ClassVar[tuple[int, ...]] = (820,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

    def __call__(self, data: ArrayLike) -> NDArray:
        rng = np.random.default_rng(self.random_seed)

        # Check if there are negative values in the data.
        values_negative = data < 0
        if values_negative.any():
            logger.warning(
                "Effect ShotNoise: %d negative pixels", values_negative.sum())
        data[values_negative] = 0

        # Apply a Poisson distribution to the low values.
        values_low = data < 1e7
        data[values_low] = rng.poisson(data[values_low])

        # Apply a normal distribution to the high values.
        values_high = ~values_low
        data[values_high] = rng.normal(data[values_high], np.sqrt(data[values_high]))

        return data

    @property
    def random_seed(self) -> int:
        return from_currsys(self.meta["random_seed"], self.cmds)

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            return det

        self.meta["random_seed"] = from_currsys(
            self.meta["random_seed"], self.cmds)

        # numpy has a problem with generating Poisson distributions above
        # certain values. E.g. on linux, numpy.random.poisson(1e20) raises
        #   ValueError: lam value too large
        # The value might be smaller on other (operating) systems.
        #
        # The poisson and normal distribution are basically the same
        # above ~100 counts:
        #   poisson(x) ~= normal(mu=x, sigma=x**0.5)
        #
        # Therefore a limit of 1e7 is used, above which the Poisson
        # distribution is approximated with a normal distribution.
        #
        # Also, the normal distribution takes only 60% as long as the
        # poisson distribution for large arrays.
        #
        # Special values should be handled with care:
        # - Negative values are mapped to 0; there cannot be negative flux.
        # - numpy.nan are implicitly passed through the normal distribution;
        #   because the Poisson distribution cannot handle them.

        new_imagehdu = fits.ImageHDU(
            data=self(det._hdu.data),
            header=det._hdu.header,
        )

        det._hdu = new_imagehdu
        return det

    def plot(self, det):
        """Plot effect image."""
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        """Plot effect histogram."""
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


class DarkCurrent(Effect):
    """
    required: dit, ndit, value
    """

    required_keys = {"value", "dit", "ndit"}
    z_order: ClassVar[tuple[int, ...]] = (830,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def __call__(
        self,
        data: np.ndarray,
        dark_level: float,
        dit: float,
        ndit: int,
    ) -> np.ndarray:
        return data + dark_level * dit * ndit

    @property
    def random_seed(self) -> int:
        return from_currsys(self.meta["random_seed"], self.cmds)

    def _get_dark_level(self, det_meta) -> float:
        dark_level = float(from_currsys(self.meta["value"], self.cmds))
        if isinstance(dark_level, Real):
            return dark_level
        if isinstance(dark_level, Mapping):
            dtcr_id = det_meta[real_colname("id", det_meta)]
            return from_currsys(dark_level[dtcr_id], self.cmds)
        raise ValueError(
            f"<{self.__class__.__name__}>.meta['value'] must be either "
            f"dict-like or scalar number, but is {dark_level}."
        )

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, Detector):
            return obj

        # Dark level needs detector meta so can't go into __call__()
        dark_level = self._get_dark_level(obj.meta)
        dit = from_currsys(self.meta["dit"], self.cmds)
        ndit = from_currsys(self.meta["ndit"], self.cmds)

        obj.data = self(obj.data, dark_level, dit, ndit)
        return obj

    def plot(self, det, **kwargs):
        """Plot effect."""
        dit = from_currsys(self.meta["dit"], self.cmds)
        ndit = from_currsys(self.meta["ndit"], self.cmds)
        total_time = dit * ndit
        times = np.linspace(0, 2*total_time, 10)
        dtcr = self.apply_to(det)
        dark_level = dtcr.data[0, 0] / total_time  # just read one pixel
        levels = dark_level * times
        fig, ax = figure_factory()
        ax.plot(times, levels, **kwargs)
        ax.set_xlabel("time")
        ax.set_ylabel("dark level")
