import pytest
from pytest import approx
import numpy as np

from scopesim.effects.electronic.noise import PoorMansHxRGReadoutNoise, _make_ron_frame
from scopesim.tests.mocks.py_objects.detector_objects import _basic_detector


class TestInit:
    def test_initialises_with_proper_keywords(self):
        isinstance(PoorMansHxRGReadoutNoise(noise_std=13, n_channels=64, ndit=1),
                   PoorMansHxRGReadoutNoise)

    def test_throws_error_without_needed_keywords(self):
        with pytest.raises(ValueError):
            PoorMansHxRGReadoutNoise()


class TestApplyTo:
    @pytest.mark.parametrize("noise_std", [13, 30, 200])
    def test_creates_noise_std_as_needed(self, noise_std):
        dtcr = _basic_detector(width=256)
        ron = PoorMansHxRGReadoutNoise(noise_std=noise_std, n_channels=64, ndit=1)
        dtcr = ron.apply_to(dtcr)

        assert np.std(dtcr._hdu.data) == approx(noise_std, rel=0.05)

    @pytest.mark.parametrize("noise_std", [1, 9, 100])
    @pytest.mark.parametrize("ndit", [1, 9, 100])
    def test_noise_reduces_with_square_root_of_ndit(self, ndit, noise_std):
        dtcr = _basic_detector(width=256)
        ron = PoorMansHxRGReadoutNoise(noise_std=noise_std, n_channels=64, ndit=ndit)
        dtcr = ron.apply_to(dtcr)
        noise_real = np.std(dtcr._hdu.data)

        assert noise_real == approx(noise_std * ndit**0.5, rel=0.1)


class TestMakeRonFrame:
    @pytest.mark.parametrize("n", (1, 4, 9, 25))
    def test_stdev_increase_with_square_of_n(self, n):
        frames = np.array([_make_ron_frame((256, 256), 10, 2, 0.1, 0.2, 0.3, 0.4)
                           for _ in range(n)])
        assert np.std(np.sum(frames, axis=0)) == approx(10*n**0.5, rel=0.3)

    @pytest.mark.parametrize("shape", [(3, 7), (7, 3)])
    def test_makes_frame_sizes_for_non_integer_n_channels(self, shape):
        n_channels = 2
        frame = _make_ron_frame(shape, 5, n_channels, 0.25, 0.25, 0.25, 0.25)
        assert frame.shape > (0, 0)
