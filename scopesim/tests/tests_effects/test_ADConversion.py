# -*- coding: utf-8 -*-
"""Contains tests for ADConversion class."""

import pytest
from unittest.mock import patch

import numpy as np

from scopesim.commands import UserCommands
from scopesim.detector import Detector, DetectorManager
from scopesim.effects.electronic import ADConversion
from scopesim.effects.detector_list import DetectorList

from scopesim.tests.mocks.py_objects.header_objects import _implane_header


@pytest.fixture(scope="class")
def patch_mock_path_micado(mock_path_micado):
    """Set the search path to the test mocks"""
    with patch("scopesim.rc.__search_path__", [mock_path_micado]):
        yield


@pytest.fixture(name="mock_detector", scope="function")
def fixture_mock_detector():
    """Instantiate a Detector object without data"""
    det = Detector(_implane_header())
    return det


@pytest.fixture(name="detector_with_data", scope="function")
def fixture_detector_with_data(mock_detector):
    """Instantiate a Detector with some data"""
    det = mock_detector
    width = det._hdu.data.shape[1]
    det._hdu.data[:] = 1.2
    det._hdu.data[:, width//2] = 1.99
    return det


# pylint: disable=missing-class-docstring,missing-function-docstring
class TestInit:
    def test_initialised_with_nothing(self):
        adconverter = ADConversion()
        assert isinstance(adconverter, ADConversion)

    def test_initialised_with_something(self):
        adconverter = ADConversion(gain=2.3, dtype="uint32")
        assert adconverter.meta['gain'] == 2.3
        assert adconverter.meta['dtype'] == "uint32"


class TestApplyTo:
    def test_converts_pixels_to_integer_with_default(self, detector_with_data):
        det = detector_with_data
        assert det.data.sum() > det.data.size * 1.2
        adconverter = ADConversion()
        adconverter.cmds = {"!DET.gain": 1, "!OBS.ndit": 1}
        det = adconverter.apply_to(det)
        assert det.data.sum() == det.data.size
        assert det.data.dtype == np.uint16

    def test_converts_pixels_to_integer_with_custom(self, detector_with_data):
        det = detector_with_data
        assert det.data.sum() > det.data.size * 1.2
        adconverter = ADConversion(dtype="int16")
        adconverter.cmds = {"!DET.gain": 1, "!OBS.ndit": 1}
        det = adconverter.apply_to(det)
        assert det.data.sum() == det.data.size
        assert det.data.dtype == np.int16

    @pytest.mark.parametrize("gain", [1.0, 2.7, 0.8])
    def test_applies_gain(self, mock_detector, gain):
        det = mock_detector
        data = np.random.randn(100, 100) * 100. + 1000.
        det._hdu.data = 1. * data
        adconverter = ADConversion(gain=gain)
        adconverter.cmds = {"!DET.gain": gain, "!OBS.ndit": 1}
        adconverter.apply_to(det)
        assert np.all((data/gain).astype(int) == det._hdu.data)

    @pytest.mark.usefixtures("patch_mock_path_micado")
    def test_applies_gain_list_to_detector_list(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0)
        detmgr = DetectorManager(det_list)
        for det in detmgr._detectors:
            det.hdu.data = 100. * np.ones_like(det.hdu.data)
        adconverter = ADConversion()
        adconverter.cmds = UserCommands(yamls=["multidetector.yaml"],
                                        properties={"!OBS.ndit": 1})
        for i, det in enumerate(detmgr._detectors):
            oldval = det.data.mean()
            newdet = adconverter.apply_to(det)
            assert newdet.data.mean() == int(oldval / (i + 1))

    @pytest.mark.parametrize("newtype, minval, maxval", [(np.uint16, 0, 65535),
                                                         (np.uint32, 0, 4294967295),
                                                         (np.int16, -32768, 32767)])
    def test_caps_data_overflow(self, mock_detector, newtype, maxval, minval):
        det = mock_detector
        medval = 100.
        gain = 2.0
        data = medval * np.ones((100, 100), dtype=np.float64)
        data[50, 50] = 1.e12
        data[60, 60] = -100000.
        det._hdu.data = data
        adconverter = ADConversion(gain=gain, dtype=newtype)
        adconverter.cmds = {"!DET.gain": gain, "!OBS.ndit": 1}
        adconverter.apply_to(det)
        assert(np.median(det._hdu.data) == newtype(medval/gain))
        assert(np.max(det._hdu.data) == maxval)
        assert(np.min(det._hdu.data) == minval)
