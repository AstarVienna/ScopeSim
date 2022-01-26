"""
Tests for Effect DetectorModePropertiesSetter
"""
import pytest
import yaml

from scopesim import rc
from scopesim import UserCommands
from scopesim.base_classes import DetectorBase
from scopesim.optics.image_plane import ImagePlane
from scopesim.effects.electronic import DetectorModePropertiesSetter
from scopesim.utils import from_currsys

from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square


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


@pytest.fixture(name="imageplane", scope="class")
def fixture_imageplane():
    """Instantiate an ImagePlane object"""
    implane = ImagePlane(_image_hdu_square().header)
    return implane


class TestInit:
    def test_initialises_with_correct_values(self):
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        assert eff.mode_properties["fast"]["!DET.mindit"] == 0.04

    def test_throws_error_when_no_dict_is_passed(self):
        with pytest.raises(ValueError):
            DetectorModePropertiesSetter()


class TestApplyTo:
    def test_currsys_updated_with_mode_specific_values(self,
                                                       imageplane):
        rc.__currsys__["!OBS.detector_readout_mode"] = "fast"
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        eff.apply_to(imageplane)

        key_name = "!DET.mindit"
        assert from_currsys(key_name) == eff.mode_properties["fast"][key_name]

    def test_understands_detector_mode_auto(self, imageplane):
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        rc.__currsys__["!OBS.auto_exposure.fill_frac"] = 0.5
        eff.apply_to(imageplane, detector_readout_mode="auto")
        assert rc.__currsys__["!OBS.detector_readout_mode"] == "slow"

    def test_throws_error_for_unknown_detector_mode(self):
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        rc.__currsys__["!OBS.detector_readout_mode"] = "notthere"
        with pytest.raises(KeyError):
            eff.apply_to(DetectorBase())

    def test_returns_object(self):
        rc.__currsys__["!OBS.detector_readout_mode"] = "fast"
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        obj = eff.apply_to(DetectorBase())

        assert isinstance(obj, DetectorBase)


class TestListModes:
    def test_it_prints(self):
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        print("\n", eff.list_modes())
        assert len(eff.list_modes()) > 0

class TestSelectModes:
    @pytest.mark.parametrize("value, mode",
                             [(10, "slow"),
                              (2e5, "middle"),
                              (1e6, "fast"),
                              (5e6, "fast")])
    def test_selects_correct_mode(self, imageplane, value, mode):
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        rc.__currsys__["!OBS.auto_exposure.fill_frac"] = 0.5
        imageplane.hdu.data += value
        assert eff.select_mode(imageplane) == mode
