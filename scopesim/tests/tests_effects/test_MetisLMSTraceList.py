"""Tests for MetisLMSSpectralTraceList"""

from unittest.mock import patch
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.interpolate import RectBivariateSpline
from astropy.io import fits
from scopesim.effects.metis_lms_trace_list import (
    _bilinear_interpolate,
    predisperser_angle,
)


# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

@pytest.fixture(scope="class")
def patch_mock_path_metis(mock_dir):
    metis_dir = mock_dir / "METIS_LMS"
    with patch("scopesim.rc.__search_path__", [metis_dir]):
        yield

class TestDetectorLayoutCache:
    def test_layout_file_is_read_only_once(self, mock_dir, monkeypatch):
        from scopesim.effects import metis_lms_trace_list as mlt

        mlt._read_detector_layout.cache_clear()
        calls = []
        real_read = mlt.ioascii.read

        def counting_read(*args, **kwargs):
            calls.append(args)
            return real_read(*args, **kwargs)

        monkeypatch.setattr(mlt.ioascii, "read", counting_read)

        path = str(mock_dir / "files" / "LIST_detector_layout.dat")
        first = mlt._read_detector_layout(path)
        second = mlt._read_detector_layout(path)

        assert len(calls) == 1
        assert first is second
        mlt._read_detector_layout.cache_clear()


class TestBilinearInterpolate:
    def test_matches_linear_spline(self):
        rng = np.random.default_rng(1234)
        cube = rng.random((5, 7, 11))
        y = rng.uniform(0, 6, (3, 4))
        x = rng.uniform(0, 10, (3, 4))

        expected = np.array([
            RectBivariateSpline(np.arange(7), np.arange(11), plane,
                                kx=1, ky=1)(y, x, grid=False)
            for plane in cube
        ])

        result = _bilinear_interpolate(cube, y, x, chunk_size=2)
        assert_allclose(result, expected, rtol=1e-14, atol=1e-14)

    def test_clips_coordinates_to_grid_boundary(self):
        cube = np.arange(24, dtype=float).reshape(2, 3, 4)
        y = np.array([[-1., 1., 3.]])
        x = np.array([[2., -1., 5.]])

        expected = np.array([
            RectBivariateSpline(np.arange(3), np.arange(4), plane,
                                kx=1, ky=1)(y, x, grid=False)
            for plane in cube
        ])

        assert_allclose(_bilinear_interpolate(cube, y, x), expected)



@pytest.mark.usefixtures("patch_mock_path_metis")
class TestPredisperserAngle:
    @pytest.mark.parametrize("coeffs,wavelen,expected",
                             [([0.], 4.5, 0),
                              ([0., 1.], 4.5, 4.5),
                              ([0., 0, 0, 1.23], 3.1, 36.64293),
                              ([-6.0585, 9.1657, -2.7017, 0.3825, -0.0205],
                               4., 6.6091)])
    def test_computes_correctly(self, coeffs, wavelen, expected):
        computed = predisperser_angle(wavelen, coeffs)
        assert_allclose(computed, expected)

    @pytest.mark.parametrize("wavelen,expected",
                             [(2.75, 5.497950),
                              (3.54, 6.280527),
                              (4.01, 6.615733),
                              (4.85, 7.138770)])
    def test_takes_coefficients_from_fits_table(self, mock_dir, wavelen,
                                                expected):
        with fits.open(mock_dir / "METIS_LMS/TRACE_LMS.fits") as hdul:
            coeff_hdu = hdul['Predisperser']
            computed = predisperser_angle(wavelen,
                                          coeff_hdu.data['coefficients'])
            assert_allclose(computed, expected, atol=1e-5)
