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
- DarkCurrent - add dark current
- LinearityCurve - apply detector (non-)linearity and saturation
- ReferencePixelBorder
- BinnedImage
- UnequalBinnedImage
- Bias - adds constant bias level to readout
"""

from ...utils import get_logger
logger = get_logger(__name__)

from .electrons import LinearityCurve, ADConversion
from .noise import (Bias, PoorMansHxRGReadoutNoise, BasicReadoutNoise,
                    ShotNoise, DarkCurrent)
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
