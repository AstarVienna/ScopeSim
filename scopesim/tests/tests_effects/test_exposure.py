"""Tests for Effect ExposureOutput"""

import pytest

import numpy as np

from scopesim import UserCommands
from scopesim.optics.image_plane import ImagePlane
from scopesim.detector import Detector
from scopesim.effects.electronic import ExposureOutput, ExposureIntegration

from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

def _patched_cmds(exptime=1, dit=None, ndit=None):
    return UserCommands(properties={"!OBS.exptime": exptime,
                                    "!OBS.dit": dit,
                                    "!OBS.ndit": ndit})

@pytest.fixture(name="imageplane", scope="class")
def fixture_imageplane():
    """Instantiate an ImagePlane object"""
    implane = ImagePlane(_image_hdu_square().header)
    implane.hdu.data += 1.e5
    return implane

@pytest.fixture(name="exposureintegration", scope="function")
def fixture_exposureintegration():
    """Instantiate an ExposureIntegration object"""
    return ExposureIntegration(dit=1., ndit=4)

@pytest.fixture(name="exposureoutput", scope="function")
def fixture_exposureoutput():
    """Instantiate an ExposureOutput object"""
    return ExposureOutput(mode="average", dit=1, ndit=4)

@pytest.fixture(name="detector", scope="function")
def fixture_detector():
    det = Detector(_image_hdu_square().header)
    det._hdu.data[:] = 1.e5
    return det

class TestExposureIntegration:
    def test_initialises_correctly(self, exposureintegration):
        assert isinstance(exposureintegration, ExposureIntegration)

    def test_integrates_correctly(self, exposureintegration, detector):
        orig = 1. * detector._hdu.data
        assert isinstance(exposureintegration.apply_to(detector), Detector)
        assert np.allclose(detector._hdu.data, 4 * orig)

    def test_fails_without_dit(self):
        with pytest.raises(ValueError):
            ExposureIntegration(ndit=5)

    def test_fails_without_ndit(self):
        with pytest.raises(ValueError):
            ExposureIntegration(dit=3.24)


class TestExposureOutput:
    def test_initialises_correctly(self, exposureoutput):
        assert isinstance(exposureoutput, ExposureOutput)

    def test_fails_with_unknown_mode(self):
        with pytest.raises(ValueError):
            ExposureOutput(mode="something", dit=1, ndit=4)

    def test_fails_without_dit_and_ndit(self):
        with pytest.raises(ValueError):
            ExposureOutput(mode="sum")

    def test_works_only_on_detector_base(self, exposureoutput, imageplane):
        assert exposureoutput.apply_to(imageplane) is imageplane

    def test_can_set_to_new_mode(self, exposureoutput):
        assert exposureoutput.current_mode == "average"
        exposureoutput.set_mode("sum")
        assert exposureoutput.current_mode == "sum"
        assert exposureoutput.meta["current_mode"] == "sum"

    def test_cannot_set_to_unknown_mode(self, exposureoutput):
        old_mode = exposureoutput.current_mode
        exposureoutput.set_mode("something")
        assert exposureoutput.current_mode == old_mode

    @pytest.mark.parametrize("dit, ndit",
                             [(1., 1),
                              (2., 5),
                              (3, 36)])
    def test_applies_average(self, dit, ndit, detector):
        det_mean = detector._hdu.data.mean()
        exposureoutput = ExposureOutput("average", dit=dit, ndit=ndit)
        result = exposureoutput.apply_to(detector)
        assert np.isclose(result._hdu.data.mean(), det_mean / ndit)

    @pytest.mark.parametrize("dit, ndit",
                             [(1., 1),
                              (2., 5),
                              (3, 36)])
    def test_applies_sum(self, dit, ndit, detector):
        det_mean = detector._hdu.data.mean()
        exposureoutput = ExposureOutput("sum", dit=dit, ndit=ndit)
        result = exposureoutput.apply_to(detector)
        assert np.isclose(result._hdu.data.mean(), det_mean)
