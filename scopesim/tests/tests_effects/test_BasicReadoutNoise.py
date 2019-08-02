import pytest
from pytest import approx
import numpy as np

from scopesim.effects import BasicReadoutNoise
from scopesim.tests.mocks.py_objects.detector_objects import _basic_detector


class TestInit:
    def test_initialises_with_proper_keywords(self):
        isinstance(BasicReadoutNoise(noise_std=13, n_channels=64, ndit=1),
                   BasicReadoutNoise)

    def test_throws_error_without_needed_keywords(self):
        with pytest.raises(ValueError):
            BasicReadoutNoise()


class TestApplyTo:
    @pytest.mark.parametrize("noise_std", [13, 30, 200])
    def test_creates_noise_std_as_needed(self, noise_std):
        dtcr = _basic_detector(width=256)
        ron = BasicReadoutNoise(noise_std=noise_std, n_channels=64, ndit=1)
        dtcr = ron.apply_to(dtcr)

        assert np.std(dtcr.image_hdu.data) == approx(noise_std, rel=0.03)

    @pytest.mark.parametrize("noise_std", [1, 9, 100])
    @pytest.mark.parametrize("ndit", [1, 9, 100])
    def test_noise_reduces_with_square_root_of_ndit(self, ndit, noise_std):
        dtcr = _basic_detector(width=256)
        ron = BasicReadoutNoise(noise_std=noise_std, n_channels=64, ndit=1)
        dtcr = ron.apply_to(dtcr)

        noise_real = np.std(dtcr.image_hdu.data)
        print(noise_real, noise_std * ndit**-0.5)

        assert noise_real == approx(noise_std * ndit**0.5, rel=0.05)
