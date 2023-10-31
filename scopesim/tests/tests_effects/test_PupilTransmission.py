import pytest
from unittest.mock import patch

from astropy import units as u

from synphot.units import PHOTLAM

from scopesim.effects.ter_curves import TERCurve
from scopesim.effects.ter_curves import PupilTransmission

# pylint: disable=no-self-use, missing-class-docstring
# pylint: disable=missing-function-docstring

THROUGHPUT = 0.764      # random value


@pytest.fixture(name="pupilmask", scope="function")
def fixture_pupilmask():
    """Instantiate a PupilTransmission object"""
    # Ensure same values no matter the currsysS
    patched = {"!SIM.spectral.wave_min": 0.3,
               "!SIM.spectral.wave_max": 20}
    with patch.dict("scopesim.rc.__currsys__", patched):
        yield PupilTransmission(transmission=THROUGHPUT)


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
