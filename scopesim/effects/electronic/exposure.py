# -*- coding: utf-8 -*-
"""Exposure actions."""

import numpy as np

from .. import Effect
from ...base_classes import DetectorBase, ImagePlaneBase
from ...utils import from_currsys, check_keys
from . import logger


class AutoExposure(Effect):
    """
    Determine DIT and NDIT automatically from ImagePlane.

    DIT is determined such that the maximum value in the incident photon flux
    (including astronomical source, sky and thermal backgrounds) fills
    the full well of the detector (``!DET.full_well``) to a given fraction
    (``!OBS.autoexposure.fill_frac``). NDIT is determined such that
    ``DIT`` * ``NDIT`` results in the requested exposure time.

    The requested exposure time is taken from ``!OBS.exptime``.

    The effects sets the parameters `!OBS.dit` and `!OBS.ndit`.

    Examples
    --------
    The parameters `!OBS.exptime`, `!DET.full_well` and `!DET.mindit` should
    be defined as properties in the respective subsections.
    ::

       name: auto_exposure
       description: automatic determination of DIT and NDIT
       class: AutoExposure
       include: True
       kwargs:
           fill_frac: "!OBS.auto_exposure.fill_frac"

    """

    required_keys = {"fill_frac", "full_well", "mindit"}

    def __init__(self, **kwargs):
        """
        The effect is the first detector effect, hence essentially operates
        on the `ImagePlane`, mapped to the detector array.
        """
        super().__init__(**kwargs)
        params = {"z_order": [902]}
        self.meta.update(params)
        self.meta.update(kwargs)
        if self.cmds is None:
            logger.error("No cmds present, using default.")
            from scopesim import UserCommands
            self.cmds = UserCommands()

        check_keys(self.meta, self.required_keys, action="error")

    def _dit_above_mindit(self, dit: float) -> bool:
        mindit = from_currsys(self.meta["mindit"], self.cmds)
        if dit < mindit:
            logger.warning("DIT = %.3f s < MINDIT = %.3f s", dit, mindit)
            return False
        return True

    def estimate_dit_ndit(
            self,
            exptime: float,
            image_plane_max: float,
            **kwargs
    ) -> tuple[float, int]:
        """
        Automatically determine DIT and NDIT from exposure time.

        Parameters
        ----------
        exptime : float
            Exposure time in seconds.
        image_plane_max : float
            Maximum pixel value from image plane, used to avoid saturation.

        Returns
        -------
        dit : float
            Detector Integration Time.
        ndit : int
            Number of Integrations.

        """
        # TODO: Remove this silly stuff once currsys works properly...
        full_well = kwargs.get(
            "full_well",
            from_currsys(self.meta["full_well"], self.cmds)
        )
        fill_frac = kwargs.get(
            "fill_frac",
            from_currsys(self.meta["fill_frac"], self.cmds)
        )

        dit_nosat = fill_frac * full_well / image_plane_max
        logger.info("Required DIT without saturation: %.3f s", dit_nosat)

        # np.ceil so that dit is at most what is required for fill_frac
        ndit = np.ceil(exptime / dit_nosat).astype(int)
        dit = exptime / ndit

        # Note: If the DIT required to avoid saturation is less than MINDIT,
        #       the observation is only possible with likely saturation...
        if not self._dit_above_mindit(dit):
            dit = from_currsys(self.meta["mindit"], self.cmds)
            # NDIT changed so that exptime is not exceeded (hence floor div)
            ndit = max(exptime // dit, 1)
            logger.warning("The detector will likely be saturated!")

        return dit, ndit

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, (ImagePlaneBase, DetectorBase)):
            return obj

        exptime = kwargs.pop("exptime",
                             from_currsys("!OBS.exptime", self.cmds))
        mindit = from_currsys(self.meta["mindit"], self.cmds)

        # TODO: Remove this silly try-except once currsys works properly...
        try:
            dit = kwargs.pop("dit", from_currsys("!OBS.dit", self.cmds))
        except (KeyError, ValueError):
            dit = None
        try:
            ndit = kwargs.pop("ndit", from_currsys("!OBS.ndit", self.cmds))
        except (KeyError, ValueError):
            ndit = None

        if dit and ndit:
            # Both DIT and NDIT are supplied (not None and non-zero), so just
            # use those regardless.
            self.cmds["!OBS.autoexpset"] = False
            # Just log warning in case DIT < MINDIT, don't actually change DIT
            self._dit_above_mindit(dit)
            return obj

        # No DIT or NDIT given, need to determine from exptime
        self.cmds["!OBS.autoexpset"] = True
        if exptime is None:
            logger.warning(
                "Please provide either !OBS.exptime or !OBS.dit + !OBS.ndit")
            if mindit is not None:
                logger.info("Using MINDIT = %.3f s for exposure time.", mindit)
                exptime = mindit
            else:
                logger.warning(
                    "MINDIT not found, falling back to 1 s for exposure time.")
                exptime = 1
        else:
            logger.info("Requested exposure time: %.3f s", exptime)

        dit, ndit = self.estimate_dit_ndit(exptime, obj.data.max(), **kwargs)

        logger.info("Exposure parameters: DIT = %.3f s, NDIT = %d", dit, ndit)
        logger.info("Total exposure time: %.3f s", dit * ndit)

        # TODO: Make sure this goes up far enough in the ChainMap...
        self.cmds["!OBS.dit"] = dit
        self.cmds["!OBS.ndit"] = ndit

        return obj


class SummedExposure(Effect):
    """Simulates a summed stack of ``ndit`` exposures."""

    required_keys = {"dit", "ndit"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [860]}
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            dit = from_currsys(self.meta["dit"], self.cmds)
            ndit = from_currsys(self.meta["ndit"], self.cmds)
            # Can't do in-place because of Quantization data type conflicts.
            obj._hdu.data = obj._hdu.data * dit * ndit

        return obj