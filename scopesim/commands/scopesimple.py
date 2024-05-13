# -*- coding: utf-8 -*-
"""Simplified top-level UI."""

from astropy.io import fits

from ..source.source import Source
from ..utils import figure_factory, top_level_catch, get_logger
from .. import OpticalTrain
from . import UserCommands


logger = get_logger(__name__)


class Simulation:
    @top_level_catch
    def __init__(self, instrument: str, **kwargs) -> None:
        self._instrument = instrument
        self._init_kwargs = kwargs
        self._last_readout = None

        self._cmds = UserCommands(use_instrument=instrument, **kwargs)
        self._opt = OpticalTrain(self._cmds)

    @top_level_catch
    def __call__(self, source: Source, **kwargs) -> fits.HDUList:
        self._opt.observe(source)
        self._last_readout = self._opt.readout(**kwargs)
        return self._last_readout[0]

    def __repr__(self) -> str:
        return f"{self.__class__}({self.instrument}, **{self._init_kwargs})"

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

    def plot(self):
        """Show simulated image."""
        _, ax = figure_factory()
        ax.imshow(self.last_readout[1], origin="lower", norm="log")
        ax.set_title(f"{self.instrument} simulation")
        return ax
