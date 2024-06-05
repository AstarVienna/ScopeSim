# -*- coding: utf-8 -*-
"""Simplified top-level UI."""

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from astropy.io import fits

from ..source.source import Source
from ..utils import figure_factory, top_level_catch, get_logger, ensure_list
from ..server import check_packages
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
    download_missing : bool
        If `instrument` isn't found, download it. Default is True, but can be
        switched off (i.e. never auto-download) by setting to False.
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
    readout()
        Re-readout the detector(s). Can be used to simulate different exptimes.

    """

    @top_level_catch
    def __init__(
            self,
            instrument: str,
            mode: str | Sequence[str] | None = None,
            download_missing: bool = True,
            **kwargs
    ) -> None:
        self._init_kwargs = kwargs
        self._last_readout = None

        check_packages(instrument, download_missing)
        if mode is not None:
            # Avoid [None], which confuses UserCommands
            mode = ensure_list(mode)
        # if "properties" in kwargs:
        #     kwargs["properties"].update({"!OBS.dit": None, "!OBS.ndit": None})
        # else:
        #     kwargs["properties"] = {"!OBS.dit": None, "!OBS.ndit": None}

        # Don't save cmds as an attribute, because OpticalTrain (currently)
        # creates a copy of that, so the actual "settings" are stored there
        # anyway. If that ever changes, consider adapting this here...
        cmds = UserCommands(use_instrument=instrument,
                            set_modes=mode, **kwargs)
        self.optical_train = OpticalTrain(cmds)

    @top_level_catch
    def __call__(
            self,
            source: Source,
            filename: Path | str | None = None,
            **kwargs
    ) -> fits.HDUList:
        """
        Run the simulation (observe the `source` and read out the detector(s)).

        Parameters
        ----------
        source : Source
            The source object to be observed by the simulation.
        filename : Path | str | None, optional
            Name or path of the FITS output file. The default is None (don't
            save anything to disk).
        **kwargs :
            Any keyword arguments passed to ``OpticalTrain.readout()``.
            Commonly used keywords are `exptime`, `dit` and `ndit`.

        Returns
        -------
        result : fits.HDUList
            List of HDUs containing simulation results.

        """
        self.optical_train.observe(source)
        return self.readout(filename, **kwargs)

    @top_level_catch
    def readout(
            self,
            filename: Path | str | None = None,
            **kwargs
    ) -> fits.HDUList:
        """
        Readout the detector(s) and optionally save the resulting FITS file.

        Only available after the simulation was run (i.e. called). This is
        usually used to run the same simulation with different exposure times,
        which can be achived by passing `exptime` as a keyword argument.

        This method is also internally called whenever the simulation is run.

        Parameters
        ----------
        filename : Path | str | None, optional
            Name or path of the FITS output file. The default is None (don't
            save anything to disk).
        **kwargs :
            Any keyword arguments passed to ``OpticalTrain.readout()``.
            Commonly used keywords are `exptime`, `dit` and `ndit`.

        Returns
        -------
        result : fits.HDUList
            List of HDUs containing simulation results.

        """
        if ("auto_exposure" in self.optical_train and
                {"dit", "ndit"}.isdisjoint(kwargs)):
            # If we have AutoExposure in the optical train and no dit or ndit
            # was passed, assume the user wants to re-estimate them from the
            # given exptime.
            kwargs["dit"] = None
            kwargs["ndit"] = None
        self._last_readout = self.optical_train.readout(filename, **kwargs)
        return self.last_readout[0]

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
        """Return instrument name.

        Shortcut for `settings["!OBS.instrument"]`.
        """
        return self.settings["!OBS.instrument"]

    @property
    def mode(self) -> str:
        """Return current instrument mode(s).

        Shortcut for `settings["!OBS.modes"]`.

        Note that some instruments use multiple "parallel" modes, in which case
        the internal list of strings is converted to a single, comma-separated
        string.
        """
        return ", ".join(self.settings["!OBS.modes"])

    @property
    def settings(self) -> UserCommands:
        """Return the current settings (UserCommands object)."""
        return self.optical_train.cmds

    def _get_vminmax(self, adjust_scale: bool):
        if not adjust_scale:
            return {"vmin": None, "vmax": None}
        imgs = self.last_readout[0][1:]
        vmin = max(min(img.data.min() for img in imgs), 0)
        vmax = max(max(img.data.max() for img in imgs), 0)
        return {"vmin": vmin, "vmax": vmax}

    @staticmethod
    def _get_figax(n_imgs: int):
        fig, axs = figure_factory(nrows=int(n_imgs), ncols=int(n_imgs))
        if isinstance(axs, np.ndarray):
            axs = list(axs.flatten())
        else:
            axs = [axs]
        return fig, axs

    def plot(self, img_slice=None, adjust_scale=False, **kwargs):
        """Show simulated image, return mpl figure and axes."""
        imgs = self.last_readout[0][1:]
        if (n_imgs := np.sqrt(len(imgs))) % 1:
            logger.warning(
                "Number of output image HDUs is not a square number, this will"
                " create an uneven image plot.")

        vminmax = self._get_vminmax(adjust_scale)
        pltkwargs = {"origin": "lower", "norm": "log"} | vminmax | kwargs
        img_slice = img_slice or slice(None)

        fig, axs = self._get_figax(n_imgs)
        for ax, img in zip(axs, imgs):
            ax.imshow(img.data[img_slice], **pltkwargs)
            ax.set_aspect("equal")
            ax.label_outer()
            ax.set_xlabel("pixel")
            ax.set_ylabel("pixel")
        fig.suptitle(f"{self.instrument} simulation in {self.mode} mode")
        return fig, axs
