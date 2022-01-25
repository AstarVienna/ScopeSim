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
            "!DET.min_dit": 0.04 
    """)


class TestInit:
    def test_initialises_with_correct_values(self):
        eff = DetectorModePropertiesSetter(**kwargs_dict())

        assert eff.mode_properties["fast"]["!DET.min_dit"] == 0.04

    def test_throws_error_when_no_dict_is_passed(self):
        with pytest.raises(ValueError):
            DetectorModePropertiesSetter()

class TestApplyTo:
    def test_currsys_properties_are_updated_with_mode_specific_values(self):
        rc.__currsys__["!OBS.detector_readout_mode"] = "fast"
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        eff.apply_to(DetectorBase())

        key_name = "!DET.min_dit"
        assert from_currsys(key_name) == eff.mode_properties["fast"][key_name]

    def test_throws_error_for_unknown_detector_mode(self):
        eff = DetectorModePropertiesSetter(**kwargs_dict())
        rc.__currsys__["!OBS.detector_readout_mode"] = "slow"
        with pytest.raises(KeyError):
            eff.apply_to(DetectorBase())
