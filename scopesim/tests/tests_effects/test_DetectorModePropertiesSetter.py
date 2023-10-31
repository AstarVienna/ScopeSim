"""
Tests for Effect DetectorModePropertiesSetter
"""
import pytest
from unittest.mock import patch

import yaml

from scopesim.base_classes import DetectorBase
from scopesim.optics.image_plane import ImagePlane
from scopesim.effects.electronic import DetectorModePropertiesSetter
from scopesim.utils import from_currsys

from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square


@pytest.fixture(scope="module")
def kwargs_dict():
    return yaml.full_load("""
    mode_properties:
        fast:
            "!DET.mindit": 0.04
            "!DET.readout_noise": 70
            "!DET.full_well": !!float 1e5
        slow:
            "!DET.mindit": 1.3
            "!DET.readout_noise": 15
            "!DET.full_well": !!float 1e5
        middle:
            "!DET.mindit": 2
            "!DET.readout_noise": 35
            "!DET.full_well": !!float 1e6
    """)


@pytest.fixture(name="basic_dmps", scope="function")
def fixture_detmodepropset(kwargs_dict):
    return DetectorModePropertiesSetter(**kwargs_dict)


@pytest.fixture(name="imageplane", scope="function")
def fixture_imageplane():
    """Instantiate an ImagePlane object"""
    implane = ImagePlane(_image_hdu_square().header)
    return implane


class TestInit:
    def test_initialises_with_correct_values(self, basic_dmps):
        assert basic_dmps.mode_properties["fast"]["!DET.mindit"] == 0.04

    def test_throws_error_when_no_dict_is_passed(self):
        with pytest.raises(ValueError):
            DetectorModePropertiesSetter()


@pytest.mark.skip(
        reason="This currently fails if run after a test using OpticalTrain, "
        "because of the rc.__currsys__ type change happening there.")
class TestApplyTo:
    @patch.dict("scopesim.rc.__currsys__",
                {"!OBS.detector_readout_mode": "fast"})
    def test_currsys_updated_with_mode_specific_values(self,
                                                       imageplane,
                                                       basic_dmps):
        basic_dmps.apply_to(imageplane)
        key_name = "!DET.mindit"
        assert from_currsys(key_name) == basic_dmps.mode_properties["fast"][key_name]

    @patch.dict("scopesim.rc.__currsys__",
                {"!OBS.detector_readout_mode": None,
                 "!OBS.auto_exposure.fill_frac": 0.5})
    def test_understands_detector_mode_auto(self, imageplane, basic_dmps):
        basic_dmps.apply_to(imageplane, detector_readout_mode="auto")
        assert from_currsys("!OBS.detector_readout_mode") == "slow"

    @patch.dict("scopesim.rc.__currsys__",
                {"!OBS.detector_readout_mode": "notthere"})
    def test_throws_error_for_unknown_detector_mode(self, basic_dmps):
        with pytest.raises(KeyError):
            basic_dmps.apply_to(DetectorBase())

    @patch.dict("scopesim.rc.__currsys__",
                {"!OBS.detector_readout_mode": "fast"})
    def test_returns_object(self, basic_dmps):
        obj = basic_dmps.apply_to(DetectorBase())
        assert isinstance(obj, DetectorBase)


class TestListModes:
    # FIXME: capture print
    def test_it_prints(self, basic_dmps):
        print("\n", basic_dmps.list_modes())
        assert len(basic_dmps.list_modes()) > 0


class TestSelectModes:
    @pytest.mark.parametrize("value, mode",
                             [(10, "slow"),
                              (2e5, "middle"),
                              (1e6, "fast"),
                              (5e6, "fast")])
    def test_selects_correct_mode(self, imageplane, value, mode, basic_dmps):
        patched = {"!OBS.auto_exposure.fill_frac": 0.5}
        with patch.dict("scopesim.rc.__currsys__", patched):
            imageplane.hdu.data += value
            assert basic_dmps.select_mode(imageplane) == mode
