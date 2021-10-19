"""Definitions of pupil masks (cold stops, etc.)"""

from astropy import units as u
from synphot.units import PHOTLAM
from ..utils import from_currsys
from .ter_curves import TERCurve

class PupilTransmission(TERCurve):
    """
    Optical element with wavelength-independent throughput

    Use this class to describe a cold stop or pupil mask that is
    characterised by "grey" throughput.
    The emissivity is set to zero, assuming that the mask is cold.
    """
    def __init__(self, throughput, **kwargs):
        wave_min = from_currsys("!SIM.spectral.wave_min") * u.um
        wave_max = from_currsys("!SIM.spectral.wave_max") *u.um
        super().__init__(wavelength=[wave_min, wave_max],
                         transmission=[throughput, throughput],
                         emissivity=[0., 0.], **kwargs)
