import os
import pytest

import numpy as np
from astropy import units as u

from synphot.units import PHOTLAM

from scopesim.effects.ter_curves import TERCurve
from scopesim.effects.ter_curves import PupilTransmission

# pylint: disable=no-self-use, missing-class-docstring
# pylint: disable=missing-function-docstring

THROUGHPUT = 0.764      # random value

@pytest.fixture(name="pupilmask", scope="class")
def fixture_pupilmask():
    """Instantiate a PupilTransmission object"""
    return PupilTransmission(transmission=THROUGHPUT)

class TestPupilTransmission:
    def test_initialises_correctly(self, pupilmask):
        assert isinstance(pupilmask, PupilTransmission)
        assert isinstance(pupilmask, TERCurve)

    def test_has_correct_throughput(self, pupilmask):
        assert pupilmask.throughput(3.5 * u.um) == THROUGHPUT

    def test_has_zero_emission(self, pupilmask):
        assert pupilmask.emission(3.5 * u.um) == 0 * PHOTLAM

    def test_update_transmission_works(self, pupilmask):
        assert pupilmask.throughput(3.5 * u.um) == THROUGHPUT
        pupilmask.update_transmission(0.5)
        assert pupilmask.throughput(3.5 * u.um) == 0.5
