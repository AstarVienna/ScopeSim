# -*- coding: utf-8 -*-
"""Simplified top-level UI."""

from collections.abc import Sequence
from pathlib import Path
import numpy as np
from astropy.io import fits

from ..source.source import Source
from ..utils import figure_factory, top_level_catch, get_logger, ensure_list
from .. import OpticalTrain
from . import UserCommands


logger = get_logger(__name__)


class Simulation:
    """
    

    Parameters
    ----------
    instrument : str
        Name of the instrument (IRDB package).
    mode : str | Sequence[str] | None, optional
        Mode or list of ombined modes for the instrument. Must match a mode
        listed in the instrument's IRDB package. If omitted, the default mode
        defined in the IRDB is used.
    **kwargs :
        Any further arguments accepted by ``scopesim.UserCommands``.

    Attributes
    ----------
    last_readout
    settings

    Methods
    -------
    plot()
        Plot the simulated image. This method can only be used once the
        instance has been called with a source object.

    """

    @top_level_catch
    def __init__(
            self,
            instrument: str,
            mode: str | Sequence[str] | None = None,
            **kwargs
    ) -> None:
        self._instrument = instrument
        self._mode = mode
        self._init_kwargs = kwargs
        self._last_readout = None

        if mode is not None:
            mode = ensure_list(mode)
        self._cmds = UserCommands(use_instrument=instrument,
                                  set_modes=mode,
                                  **kwargs)
        self.optical_train = OpticalTrain(self._cmds)

    @top_level_catch
    def __call__(
            self,
            source: Source,
            filename: Path | str | None = None,
            **kwargs
    ) -> fits.HDUList:
        self.optical_train.observe(source)
        self._last_readout = self.optical_train.readout(filename, **kwargs)
        return self._last_readout[0]

    def __repr__(self) -> str:
        return (f"{self.__class__}({self.instrument}, {self.mode}, "
                f"**{self._init_kwargs})")

    def __str__(self) -> str:
        return f"ScopeSim Simulation for {self.instrument} in {self.mode} mode"

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
        """Show simulated image, return mpl figure and axes."""
        imgs = self.last_readout[0][1:]
        if (n_imgs := np.sqrt(len(imgs))) % 1:
            logger.warning(
                "Number of output image HDUs is not a square number, this will"
                " create an uneven image plot.")
        vmin = max(min(img.data.min() for img in imgs), 0)
        vmax = max(max(img.data.max() for img in imgs), 0)
        fig, axs = figure_factory(nrows=int(n_imgs), ncols=int(n_imgs))
        if isinstance(axs, np.ndarray):
            axs = axs.flatten()
        else:
            axs = [axs]
        for ax, img in zip(axs, imgs):
            ax.imshow(img.data, origin="lower",
                      norm="log", vmin=vmin, vmax=vmax)
            ax.set_aspect("equal")
            ax.label_outer()
        fig.suptitle(f"{self.instrument} simulation")
        return fig, axs
