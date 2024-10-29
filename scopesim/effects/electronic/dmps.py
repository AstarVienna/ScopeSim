# -*- coding: utf-8 -*-
"""Outdated and under consideration for removal."""

from typing import ClassVar

import numpy as np

from .. import Effect
from ...base_classes import ImagePlaneBase
from ...utils import from_currsys, check_keys, pretty_print_dict
from . import logger


class DetectorModePropertiesSetter(Effect):
    """
    Set mode specific curr_sys properties for different detector readout modes.

    A little class (``DetectorModePropertiesSetter``) that allows different
    ``"!DET"`` properties to be set on the fly.

    Parameters
    ----------
    mode_properties : dict
        A dictionary containing the DET parameters to be changed for each mode.
        See below for an example yaml entry.

    Examples
    --------
    Add the values for the different detector readout modes to all the relevant
    detector yaml files. In this case the METIS HAWAII (L, M band) and GeoSnap
    (N band) detectors: METIS_DET_IMG_LM.yaml , METIS_DET_IMG_N.yaml
    ::

        - name: lm_detector_readout_parameters
          class: DetectorModePropertiesSetter
          kwargs:
            mode_properties:
              fast:
                mindit: 0.04
                full_well: !!float 1e5
                ron: 70
              slow:
                mindit: 1.3
                full_well: !!float 1e5
                ron: 14

    Add the OBS dict entry !OBS.detector_readout_mode to the properties section
    of the mode_yamls descriptions in the default.yaml files.
    ::

        mode_yamls:
          - object: observation
            alias: OBS
            name: lss_l
            yamls:
              ...
            properties:
              ...
              detector_readout_mode: slow

    """

    required_keys = {"mode_properties"}
    z_order: ClassVar[tuple[int, ...]] = (299, 900)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

        self.mode_properties = kwargs["mode_properties"]

    def apply_to(self, obj, **kwargs):
        mode_name = kwargs.get("detector_readout_mode",
                               from_currsys("!OBS.detector_readout_mode",
                                            self.cmds))
        if isinstance(obj, ImagePlaneBase) and mode_name == "auto":
            mode_name = self.select_mode(obj, **kwargs)
            logger.info("Detector mode set to %s", mode_name)

        self.meta["detector_readout_mode"] = mode_name
        props_dict = self.mode_properties[mode_name]
        self.cmds["!OBS.detector_readout_mode"] = mode_name
        for key, value in props_dict.items():
            self.cmds[key] = value

        return obj

    def list_modes(self):
        """Return list of available detector modes."""
        return pretty_print_dict(self.mode_properties)

    def select_mode(self, obj, **kwargs):
        """Automatically select detector mode based on image plane peak value.

        Select the mode with lowest readnoise that does not saturate the
        detector. When all modes saturate, select the mode with the lowest
        saturation level (peak to full_well).
        """
        immax = np.max(obj.data)
        fillfrac = kwargs.get("fill_frac",
                              from_currsys("!OBS.auto_exposure.fill_frac",
                                           self.cmds))

        goodmodes = []
        goodron = []
        badmodes = []
        satlevel = []
        for modeid, modeprops in self.mode_properties.items():
            mindit = modeprops["!DET.mindit"]
            fullwell = modeprops["!DET.full_well"]
            if immax * mindit < fillfrac * fullwell:
                goodmodes.append(modeid)
                goodron.append(modeprops["!DET.readout_noise"])
            else:
                badmodes.append(modeid)
                satlevel.append(immax * mindit / fullwell)

        if not goodmodes:
            return badmodes[np.argmin(satlevel)]

        return goodmodes[np.argmin(goodron)]
