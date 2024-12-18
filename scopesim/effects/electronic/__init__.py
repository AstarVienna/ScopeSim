# -*- coding: utf-8 -*-
"""
Electronic detector effects - related to detector readout.

Classes:
- DetectorModePropertiesSetter - set parameters for readout mode
- AutoExposure - determine DIT and NDIT automatically
- SummedExposure - simulates a summed stack of ``ndit`` exposures
- PoorMansHxRGReadoutNoise - simple readout noise for HAWAII detectors
- BasicReadoutNoise - readout noise
- ShotNoise - realisation of Poissonian photon noise
- DarkCurrent - add dark current
- LinearityCurve - apply detector (non-)linearity and saturation
- ReferencePixelBorder
- BinnedImage
- UnequalBinnedImage
- Bias - adds constant bias level to readout
"""

from ...utils import get_logger
logger = get_logger(__name__)

from .electrons import LinearityCurve, Quantization
from .noise import (Bias, PoorMansHxRGReadoutNoise, BasicReadoutNoise,
                    ShotNoise, DarkCurrent)
from .exposure import AutoExposure, SummedExposure
from .pixels import ReferencePixelBorder, BinnedImage, UnequalBinnedImage
from .dmps import DetectorModePropertiesSetter
