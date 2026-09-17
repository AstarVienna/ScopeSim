# -*- coding: utf-8 -*-
"""Any kinds of electronic or photonic noise."""

from hashlib import sha256
from typing import ClassVar
from collections.abc import Mapping
from numbers import Real  # matches int, float and all the numpy scalars

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ...utils import from_currsys, figure_factory
from . import ElectronicEffect, Detector, logger


# TODO: Potential refactoring in multiple effects here:
#       The various __call__ now get a det_id, which is what the _get_... need
#       for their lookup. Either pass the rng to __call__ instead of det_id,
#       or call the _get_.. in there. Needs some consideration of outside
#       callers, would they rather pass which one? And how about det_id = None?
#       Also could refactor the _get_... (also elsewhere) to a abstract value
#       plus the concrete key.


class Bias(ElectronicEffect):
    """Adds a constant bias level to readout."""

    required_keys = {"bias"}
    z_order: ClassVar[tuple[int, ...]] = (855,)

    def __call__(self, data: ArrayLike) -> NDArray:
        return data + self.bias_level

    @property
    def bias_level(self) -> float:
        return from_currsys(self.meta["bias"], self.cmds)


class DarkCurrent(ElectronicEffect):
    """
    required: dit, ndit, value
    """

    required_keys = {"value", "dit", "ndit"}
    z_order: ClassVar[tuple[int, ...]] = (830,)

    def __call__(
        self,
        data: ArrayLike,
        dark_level: float,
        dit: float,
        ndit: int,
    ) -> NDArray:
        return data + dark_level * dit * ndit

    def _get_dark_level(self, det_id: int) -> float:
        dark_level = float(from_currsys(self.meta["value"], self.cmds))
        if isinstance(dark_level, Real):
            return dark_level
        if isinstance(dark_level, Mapping):
            return from_currsys(dark_level[det_id], self.cmds)
        raise TypeError(
            f"<{self.__class__.__name__}>.meta['value'] must be either "
            f"dict-like or scalar number, but is {dark_level}."
        )

    def _apply_to_det(self, det: Detector) -> None:
        logger.debug("Apply %s to %s", self.display_name, det)

        # Dark level needs detector meta so can't go into __call__()
        dark_level = self._get_dark_level(det.det_id)
        dit = from_currsys(self.meta["dit"], self.cmds)
        ndit = from_currsys(self.meta["ndit"], self.cmds)

        det.data = self(det.data, dark_level, dit, ndit)

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
        return ax


# TODO: Add tests for seed resolving, start with code from #975, but for
#       basic_instrument, and work from there.
class RandomEffect(ElectronicEffect):
    """Mixin class for Effects that need random seeds."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["random_seed"] = "!SIM.random.seed"  # Default

    @property
    def root_seed(self) -> int:
        """Resolve root seed from cmds, or use value from effect kwargs."""
        return from_currsys(self.meta["random_seed"], self.cmds)

    @classmethod
    def cls_seed(cls) -> int:
        """Generate reproducible seed from hash of class name."""
        digest = sha256(cls.__name__.encode("utf-8")).digest()
        return int.from_bytes(digest[:16], "little")

    @property
    def readout_seed(self) -> int:
        """Readout ID, or 0 if not found in cmds."""
        return self.cmds.get("!OBS.roid", 0) if self.cmds is not None else 0

    @property
    def random_seed(self) -> int | None:
        """Composite random seed, or None."""
        seed = [self.readout_seed, self.cls_seed(), self.root_seed]
        # TODO: Consider removing this if root seed is resolved upstream!
        if None in seed:
            return None
        return seed

    def create_rng(self, det_id: int | None = None) -> np.random.Generator:
        """
        Instantiate np.random.Generator using composite seed.

        Parameters
        ----------
        det_id : int | None, optional
            Detector ID. If None (the default), detector ID is not included in
            the composite seed, meaning results will look identical on all
            detectors passed to the effect in the same readout.

        Returns
        -------
        np.random.Generator
            New Generator instance.

        """
        if det_id is not None:
            seed = [det_id, *self.random_seed]
        else:
            seed = self.random_seed
        return np.random.default_rng(seed)

    def plot(self, det):
        """Plot effect image."""
        detector = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(detector.data)
        return ax

    def plot_hist(self, det, **kwargs):
        """Plot effect histogram."""
        detector = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(detector.data.flatten())
        return ax


class BasicReadoutNoise(RandomEffect):
    """Readout noise computed as: ron * sqrt(NDIT)."""

    required_keys = {"noise_std", "ndit"}
    z_order: ClassVar[tuple[int, ...]] = (811,)

    def __call__(
        self,
        data: ArrayLike,
        ndit: int = 1,
        det_id: int | None = None,
    ) -> NDArray:
        rng = self.create_rng(det_id)
        return data + self._create_noise_frame(data.shape, rng, ndit)

    @property
    def noise_std(self) -> float:
        return from_currsys(self.meta["noise_std"], self.cmds)

    def _create_noise_frame(
        self,
        shape: tuple[int, ...],
        rng: np.random.Generator,
        ndit: int,
    ) -> NDArray:
        scale = self.noise_std * np.sqrt(float(ndit))
        return rng.normal(loc=0, scale=scale, size=shape)

    def _apply_to_det(self, det: Detector) -> None:
        logger.debug("Apply %s to %s", self.display_name, det)

        ndit = from_currsys(self.meta["ndit"], self.cmds)
        det.data = self(det.data, ndit)


class PoorMansHxRGReadoutNoise(BasicReadoutNoise):
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
        }
        self.meta.update(params)
        self.meta.update(kwargs)

    @property
    def n_channels(self) -> int:
        return from_currsys(self.meta["n_channels"], self.cmds)

    def _create_noise_frame(
        self,
        shape: tuple[int, ...],
        rng: np.random.Generator,
        ndit: int,
    ) -> NDArray:
        ron_frame = self._make_ron_frame(rng, shape)
        stacked_ron_frame = np.zeros_like(ron_frame)

        for i in range(ndit):
            stacked_ron_frame += np.roll(
                ron_frame,
                rng.integers((0, 0), ron_frame.shape),
                axis=(0, 1),
            )

        return stacked_ron_frame

    def _make_ron_frame(self, rng, shape: tuple[int, ...]) -> NDArray:
        self.meta = from_currsys(self.meta, self.cmds)  # TODO: Is this needed?
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

        return (pixel + line) + channel[:shape[0], None]

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


class PixelResponseNonUniformity(RandomEffect):
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

    include:  "!DET.include_prnu"

    Example
    -------
    ::

       - name: prnu
         description: Pixel response non-uniformity
         class: PixelResponseNonUniformity
         kwargs:
           prnu_std: 0.001
           include: "!DET.include_prnu"

    """

    required_keys: ClassVar[set] = {"prnu_std"}
    z_order: ClassVar[tuple[int, ...]] = (805,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._gain_maps = {}  # keyed by det_id

    def __call__(
        self,
        data: ArrayLike,
        prnu_std: float,
        det_id: int | None = None,
    ) -> NDArray:
        if det_id not in self._gain_maps:
            # TODO: Should not use roid here???
            rng = self.create_rng(det_id)
            self._gain_maps[det_id] = rng.normal(
                loc=1.0, scale=prnu_std, size=data.shape,
            )

        return data * self._gain_maps[det_id]

    @property
    def random_seed(self) -> int | None:
        """Composite random seed, or None.

        Override from base class to remove ROID, because PRNU should be
        identical in one observation.
        """
        seed = [self.cls_seed(), self.root_seed]
        # TODO: Consider removing this if root seed is resolved upstream!
        if None in seed:
            return None
        return seed

    def _get_prnu_std(self, det_id: int) -> float:
        prnu_std = from_currsys(self.meta["prnu_std"], self.cmds)
        if isinstance(prnu_std, Real):
            return prnu_std
        if isinstance(prnu_std, Mapping):
            return from_currsys(prnu_std[det_id], self.cmds)
        raise TypeError(
            f"<{self.__class__.__name__}>.meta['value'] must be either "
            f"dict-like or scalar number, but is {prnu_std}."
        )

    def _apply_to_det(self, det: Detector) -> None:
        logger.debug("Apply %s to %s", self.display_name, det)

        prnu_std = self._get_prnu_std(det.det_id)
        det.data = self(det.data, prnu_std)

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
        return fig, ax


class ShotNoise(RandomEffect):
    """Poissonian photon noise.

    Notes
    -----

    Numpy has a problem with generating Poisson distributions above certain
    values. E.g. on linux, numpy.random.poisson(1e20) raises ValueError: lam
    value too large. The value might be smaller on other (operating) systems.

    The poisson and normal distribution are basically the same
    above ~100 counts:
      poisson(x) ~= normal(mu=x, sigma=x**0.5)

    Therefore a limit of 1e7 is used, above which the Poisson distribution is
    approximated with a normal distribution.

    Also, the normal distribution takes only 60% as long as the Poisson
    distribution for large arrays.

    Special values should be handled with care:
    - Negative values are mapped to 0; there cannot be negative flux.
    - numpy.nan are implicitly passed through the normal distribution;
      because the Poisson distribution cannot handle them.
    """
    z_order: ClassVar[tuple[int, ...]] = (820,)

    def __call__(self, data: ArrayLike, det_id: int | None = None) -> NDArray:
        rng = self.create_rng(det_id)

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
        data[values_high] = rng.normal(
            loc=data[values_high],
            scale=np.sqrt(data[values_high]),
        )

        return data
