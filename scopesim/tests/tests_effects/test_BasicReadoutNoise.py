import pytest
from pytest import approx
import numpy as np

from scopesim.effects import BasicReadoutNoise, make_ron_frame
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

        assert np.std(dtcr._hdu.data) == approx(noise_std, rel=0.05)

    @pytest.mark.parametrize("noise_std", [1, 9, 100])
    @pytest.mark.parametrize("ndit", [1, 9, 100])
    def test_noise_reduces_with_square_root_of_ndit(self, ndit, noise_std):
        dtcr = _basic_detector(width=256)
        ron = BasicReadoutNoise(noise_std=noise_std, n_channels=64, ndit=ndit)
        dtcr = ron.apply_to(dtcr)
        noise_real = np.std(dtcr._hdu.data)

        assert noise_real == approx(noise_std * ndit**0.5, rel=0.05)


class TestMakeRonFrame:
    @pytest.mark.parametrize("n", (1, 4, 9, 25))
    def test_stdev_increase_with_square_of_n(self, n):
        frames = np.array([make_ron_frame((256, 256), 10, 2, 0.1, 0.2, 0.3, 0.4)
                           for _ in range(n)])
        assert np.std(np.sum(frames, axis=0)) == approx(10*n**0.5, rel=0.1)
