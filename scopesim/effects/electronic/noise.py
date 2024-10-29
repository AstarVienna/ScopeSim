# -*- coding: utf-8 -*-
"""Any kinds of electronic or photonic noise."""

from typing import ClassVar

import numpy as np
from astropy.io import fits

from .. import Effect
from ...base_classes import DetectorBase
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

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            biaslevel = from_currsys(self.meta["bias"], self.cmds)
            # Can't do in-place because of Quantization data type conflicts.
            obj._hdu.data = obj._hdu.data + biaslevel

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

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"],
                                                    self.cmds)
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            self.meta = from_currsys(self.meta, self.cmds)
            ron_keys = ["noise_std", "n_channels", "channel_fraction",
                        "line_fraction", "pedestal_fraction", "read_fraction"]
            ron_kwargs = {key: self.meta[key] for key in ron_keys}
            ron_kwargs["image_shape"] = det._hdu.data.shape

            ron_frame = _make_ron_frame(**ron_kwargs)
            stacked_ron_frame = np.zeros_like(ron_frame)
            for i in range(self.meta["ndit"]):
                dx = np.random.randint(0, ron_frame.shape[1])
                dy = np.random.randint(0, ron_frame.shape[0])
                stacked_ron_frame += np.roll(ron_frame, (dy, dx), axis=(0, 1))

            # TODO: this .T is ugly. Work out where things are getting switched and remove it!
            # Can't do in-place because of Quantization data type conflicts.
            det._hdu.data = det._hdu.data + stacked_ron_frame.T

        return det

    def plot(self, det, **kwargs):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data, origin="lower")

    def plot_hist(self, det, **kwargs):
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

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            ndit = from_currsys(self.meta["ndit"], self.cmds)
            ron = from_currsys(self.meta["noise_std"], self.cmds)
            noise_std = ron * np.sqrt(float(ndit))

            random_seed = from_currsys(self.meta["random_seed"], self.cmds)
            if random_seed is not None:
                np.random.seed(random_seed)
            # Can't do in-place because of Quantization data type conflicts.
            det._hdu.data = det._hdu.data + np.random.normal(
                loc=0, scale=noise_std, size=det._hdu.data.shape)

        return det

    def plot(self, det):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


class ShotNoise(Effect):
    z_order: ClassVar[tuple[int, ...]] = (820,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

    def apply_to(self, det, **kwargs):
        if not isinstance(det, DetectorBase):
            return det

        self.meta["random_seed"] = from_currsys(self.meta["random_seed"],
                                                self.cmds)
        rng = np.random.default_rng(self.meta["random_seed"])

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

        data = det._hdu.data

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

        new_imagehdu = fits.ImageHDU(
            data=data,
            header=det._hdu.header,
        )

        det._hdu = new_imagehdu
        return det

    def plot(self, det):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
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

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            if isinstance(from_currsys(self.meta["value"], self.cmds), dict):
                dtcr_id = obj.meta[real_colname("id", obj.meta)]
                dark = from_currsys(self.meta["value"][dtcr_id], self.cmds)
            elif isinstance(from_currsys(self.meta["value"], self.cmds), float):
                dark = from_currsys(self.meta["value"], self.cmds)
            else:
                raise ValueError("<DarkCurrent>.meta['value'] must be either "
                                 f"dict or float, but is {self.meta['value']}")

            dit = from_currsys(self.meta["dit"], self.cmds)
            ndit = from_currsys(self.meta["ndit"], self.cmds)

            # Can't do in-place because of Quantization data type conflicts.
            obj._hdu.data = obj._hdu.data + dark * dit * ndit

        return obj

    def plot(self, det, **kwargs):
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


def _pseudo_random_field(scale=1, size=(1024, 1024)):
    n = 256
    image = np.zeros(size)
    batch = np.random.normal(loc=0, scale=scale, size=(2*n, 2*n))
    for y in range(0, size[1], n):
        for x in range(0, size[0], n):
            i, j = np.random.randint(n, size=2)
            dx, dy = min(size[0]-x, n), min(size[1]-y, n)
            image[x:x+dx, y:y+dy] = batch[i:i+dx, j:j+dy]

    return image


def _make_ron_frame(image_shape, noise_std, n_channels, channel_fraction,
                   line_fraction, pedestal_fraction, read_fraction):
    shape = image_shape
    w_chan = max(1, shape[0] // n_channels)

    pixel_std = noise_std * (pedestal_fraction + read_fraction)**0.5
    line_std = noise_std * line_fraction**0.5
    if shape < (1024, 1024):
        pixel = np.random.normal(loc=0, scale=pixel_std, size=shape)
        line = np.random.normal(loc=0, scale=line_std, size=shape[1])
    else:
        pixel = _pseudo_random_field(scale=pixel_std, size=shape)
        line = pixel[0]

    channel_std = noise_std * channel_fraction**0.5
    channel = np.repeat(np.random.normal(loc=0, scale=channel_std,
                                         size=n_channels), w_chan + 1, axis=0)

    ron_frame = (pixel + line).T + channel[:shape[0]]

    return ron_frame
