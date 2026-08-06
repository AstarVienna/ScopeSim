"""Unit tests for mosaic_trace_list.py"""

# pylint: disable=missing-function-docstring
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
from unittest.mock import patch
import pytest

import numpy as np
from astropy.io import fits
from scopesim.effects.mosaic_trace_list import (Transform1D,
                                                MosaicSpectralTraceList,
                                                MosaicSpectralTrace,
                                                MosaicCollapseSpectralTraces,
                                                MosaicConvertToTable)
from scopesim.optics.image_plane_utils import header_from_list_of_xy
from scopesim.detector import Detector

@pytest.fixture(name="tf1d", scope="class")
def fixture_tf1d():
    """Instantiate a Transform1D"""
    coeffs = np.array([2, -1, 1])
    return Transform1D(coeffs)

@pytest.fixture(name="quadratic", scope="class")
def fixture_quadratic():
    """Quadratic model, analytic and coeffients"""
    coeffs = np.array([1, -1, 2])

    def quadfunc(x):
        z_a = 1 - 1 * x + 2 * x**2
        return z_a

    def dquad_dx(x):
        return -1 + 4 * x

    return {'coeffs': coeffs,
            'function': quadfunc,
            'gradient': dquad_dx}

class TestTransform1D:
    """Tests for Transform1D()"""
    def test_initialises_with_coeffs(self, tf1d):
        assert isinstance(tf1d, Transform1D)

    def test_call_gives_correct_result(self, quadratic):
        x = np.random.randn()

        # coefficients and explicit function
        tf1d = Transform1D(quadratic['coeffs'])
        assert tf1d(x) == quadratic['function'](x)

    def test_gradient_gives_correct_result(self, quadratic):
        x = np.random.randn()

        tf2d = Transform1D(quadratic['coeffs'])
        tf2d_grad = tf2d.gradient()

        assert tf2d_grad(x) == quadratic['gradient'](x)

    def test_fit_gives_correct_coeffs(self):
        x = np.linspace(0, 1, 10)
        y = 1 - 0.5 * x + 2.3 * x**2 - 3 * x**3

        coeffs = np.array([1, -0.5, 2.3, -3])
        tf1d = Transform1D.fit(x, y, degree=3)

        assert tf1d.coeffs == pytest.approx(coeffs)

@pytest.fixture(scope="class")
def patch_mock_path_mosaic(mock_dir):
    mosaic_dir = mock_dir / "MOSAIC"
    with patch("scopesim.rc.__search_path__", [mosaic_dir]):
        yield

@pytest.fixture(name="tracelist", scope="class")
def fixture_tracelist():
    """Instantiate a MosaicSpectralTraceList"""
    return MosaicSpectralTraceList(filename="TRACE_MOS-LR-J.fits")

@pytest.fixture(name="collapse1d", scope="class")
def fixture_collapse1d():
    """Instantiate a CollapseSpectralTrace"""
    return MosaicCollapseSpectralTraces(filename="TRACE_MOS-LR-J.fits")

@pytest.fixture(name="convert2table", scope="class")
def fixture_convert2table():
    """Instantiate a ConvertToTable"""
    return MosaicConvertToTable(filename="TRACE_MOS-LR-J.fits")

@pytest.fixture(name="detector", scope="class")
def fixture_detector():
    """Instantiate a detector image"""
    rows = np.array([2023, 2031, 2039, 2047, 2055, 2063, 2071])
    xsize, ysize = 50, 4096
    hdr = header_from_list_of_xy([-xsize/2, xsize/2], [-ysize/2, ysize/2], 1, "D")
    dtcr = Detector(hdr)
    dtcr._hdu.data[rows, :] = 1
    return dtcr

@pytest.fixture(name="apply_collapse1d", scope="class")
def fixture_apply_collapse1d(collapse1d, detector):
    return collapse1d.apply_to(detector)

@pytest.mark.usefixtures("patch_mock_path_mosaic")
class TestSpectralTraceList:
    """Tests for MosaicSpectralTraceList"""
    def test_initialise_tracelist(self, tracelist):
        assert isinstance(tracelist, MosaicSpectralTraceList)

    def test_tracelist_loads_aperture_list(self, tracelist):
        assert len(tracelist.aplist) == 7
        assert "right" in tracelist.aplist.colnames

    def test_tracelist_creates_traces(self, tracelist):
        assert len(tracelist.spectral_traces) == 7

    def test_trace_has_transforms(self, tracelist):
        for trace in tracelist.spectral_traces.values():
            assert isinstance(trace.lam2x, Transform1D)
            assert isinstance(trace.lam2y, Transform1D)
            assert isinstance(trace.x2lam, Transform1D)


@pytest.mark.usefixtures("patch_mock_path_mosaic")
class TestOutputFormats:
    """Tests for MOSAIC MOS/mIFU output formats"""
    def test_initialise_collapse1d(self, collapse1d):
        assert isinstance(collapse1d, MosaicCollapseSpectralTraces)

    def test_collapse1d_returns_detector(self, apply_collapse1d):
        assert isinstance(apply_collapse1d, Detector)

    def test_collapse1d_detector_has_bintable(self, apply_collapse1d):
        assert isinstance(apply_collapse1d._hdu, fits.BinTableHDU)

    def test_collapse1d_gives_correct_result(self, apply_collapse1d):
        assert np.all(apply_collapse1d._hdu.data['spectrum'] == 7)

    def test_initialise_table(self, convert2table):
        assert isinstance(convert2table, MosaicConvertToTable)
