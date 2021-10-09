"""Definitions of pupil masks (cold stops, etc.)"""

from synphot.units import PHOTLAM
from .ter_curves import TERCurve

class PupilTransmission(TERCurve):
    """
    Optical element with wavelength-independent throughput

    Use this class to describe a cold stop or pupil mask that is
    characterised by "grey" throughput.
    The emission has been set to zero, assuming that the mask is cold.
    """
    def __init__(self, throughput, **kwargs):
        super().__init__(wavelength=[1., 2.],
                         transmission=[throughput, throughput],
                         emissivity=[0., 0.])
