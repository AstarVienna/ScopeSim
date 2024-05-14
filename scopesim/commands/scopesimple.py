# -*- coding: utf-8 -*-
"""Simplified top-level UI."""

import numpy as np
from astropy.io import fits

from ..source.source import Source
from ..utils import figure_factory, top_level_catch, get_logger, ensure_list
from .. import OpticalTrain
from . import UserCommands


logger = get_logger(__name__)


class Simulation:
    @top_level_catch
    def __init__(self, instrument: str, mode: str, **kwargs) -> None:
        self._instrument = instrument
        self._mode = mode
        self._init_kwargs = kwargs
        self._last_readout = None

        self._cmds = UserCommands(use_instrument=instrument,
                                  set_modes=ensure_list(mode),
                                  **kwargs)
        self.optical_train = OpticalTrain(self._cmds)

    @top_level_catch
    def __call__(self, source: Source, **kwargs) -> fits.HDUList:
        self.optical_train.observe(source)
        self._last_readout = self.optical_train.readout(**kwargs)
        return self._last_readout[0]

    def __repr__(self) -> str:
        return (f"{self.__class__}({self.instrument}, {self.mode}, "
                f"**{self._init_kwargs})")

    def __str__(self) -> str:
        return f"ScopeSim Simulation for {self.instrument}"

    def _repr_pretty_(self, printer, cycle):
        """For nice ipython display."""
        if cycle:
            printer.text("Simulation(...)")
        else:
            printer.text(str(self))

    @property
    def last_readout(self) -> fits.HDUList:
        """Return last readout HDUL, if any."""
        if self._last_readout is None:
            logger.error("No readout found, run simulation first!")
            return
        return self._last_readout

    @property
    def instrument(self) -> str:
        """Return instrument name."""
        return self._instrument

    @property
    def mode(self) -> str:
        """Return current instrument mode."""
        return self._mode

    @property
    def settings(self) -> UserCommands:
        """Return the current settings (UserCommands object)."""
        return self._cmds

    def plot(self):
        """Show simulated image."""
        imgs = self.last_readout[1:]
        if (n_imgs := np.sqrt(len(imgs))) % 1:
            logger.warning(
                "Number of output image HDUs is not a square number, this will"
                " create an uneven image plot.")
        fig, axs = figure_factory(nrows=int(n_imgs), ncols=int(n_imgs))
        for ax, img in zip(axs, imgs):
            ax.imshow(img.data, origin="lower", norm="log")
        fig.set_title(f"{self.instrument} simulation")
        return axs
