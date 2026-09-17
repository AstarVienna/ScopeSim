# -*- coding: utf-8 -*-
"""
Electronic detector effects - related to detector readout.

Classes:
- DetectorModePropertiesSetter - set parameters for readout mode
- AutoExposure - determine DIT and NDIT automatically
- ExposureIntegration - integrates flux over exposure time
- PoorMansHxRGReadoutNoise - simple readout noise for HAWAII detectors
- BasicReadoutNoise - readout noise
- ShotNoise - realisation of Poissonian photon noise
- PixelResponseNonUniformity - per-pixel gain variation (PRNU)
- DarkCurrent - add dark current
- LinearityCurve - apply detector (non-)linearity and saturation
- ReferencePixelBorder
- BinnedImage
- UnequalBinnedImage
- Bias - adds constant bias level to readout
- InterPixelCapacitance - apply IPC kernel to detector readout
"""

from ...utils import get_logger, check_keys
from ...detector import Detector
from .. import Effect

logger = get_logger(__name__)


class ElectronicEffect(Effect):
    """Base class for electronic effects.

    This will eventually replace the 800-range zorder (or parts of it).
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)
        check_keys(self.meta, self.required_keys, action="error")

    def _apply_to_det(self, det: Detector) -> None:
        """Subclasses can override if more params needed in call."""
        logger.debug("Apply %s to %s", self.display_name, det)
        det.data = self(det.data)

    def apply_to(self, obj, **kwargs):
        """See parent docstring."""
        if isinstance(obj, Detector):
            self._apply_to_det(obj)

        return obj


from .electrons import LinearityCurve, ADConversion, InterPixelCapacitance
from .noise import (Bias, PoorMansHxRGReadoutNoise, BasicReadoutNoise,
                    ShotNoise, DarkCurrent, PixelResponseNonUniformity)
from .exposure import AutoExposure, ExposureIntegration, ExposureOutput
from .pixels import ReferencePixelBorder, BinnedImage, UnequalBinnedImage
from .dmps import DetectorModePropertiesSetter


# TODO: rm this in v1.0
def Quantization(*args, **kwargs):
    raise AttributeError(
        "The `Quantization` effect was removed in v0.10.0. Please update the "
        "requested IRDB package by running `download_packages(<package_name>)`"
        "or by updating your local IRDB clone.")

def SummedExposure(*args, **kwargs):
    raise AttributeError(
        "The `SummedExposure` effect was removed in v0.10.0. Please update the "
        "requested IRDB package by running `download_packages(<package_name>)`"
        "or by updating your local IRDB clone.")
