"""Tests for MetisLMSSpectralTraceList"""

from unittest.mock import patch
import pytest
from numpy.testing import assert_allclose
from astropy.io import fits
from scopesim.effects.metis_lms_trace_list import predisperser_angle


# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

@pytest.fixture(scope="class")
def patch_mock_path_metis(mock_dir):
    metis_dir = mock_dir / "METIS_LMS"
    with patch("scopesim.rc.__search_path__", [metis_dir]):
        yield

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
