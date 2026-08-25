"""Tests for MetisLMSSpectralTraceList"""

from unittest.mock import patch
import pytest
from numpy.testing import assert_allclose
from astropy.io import fits
from scopesim.effects.metis_lms_trace_list import (predisperser_angle,
                                                   echelle_setting,
                                                   MetisLMSSpectralTrace)


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

class TestMetisLMSSpectralTraceInit:
    def test_does_not_mutate_caller_params(self, mock_dir, monkeypatch):
        """The params dict is passed by reference from the trace list's
        self.meta; the trace previously wrote its slice-specific values
        (slice, aperture_id, ...) back into it."""
        monkeypatch.setattr(MetisLMSSpectralTrace, "fov_grid",
                            lambda self: {})

        with fits.open(mock_dir / "METIS_LMS/TRACE_LMS.fits") as hdul:
            ech = echelle_setting(4.2, 18.2, hdul["WCAL"].data)
            params = {"order": ech["Ord"], "echelle": ech["Echelle"],
                      "wavelen": 4.2}
            snapshot = dict(params)

            trace = MetisLMSSpectralTrace(hdul, spslice=0, params=params)

            assert params == snapshot, \
                "constructor modified the caller's params dict"
            assert trace.meta["slice"] == 0
            assert trace.meta["aperture_id"] == 0


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
