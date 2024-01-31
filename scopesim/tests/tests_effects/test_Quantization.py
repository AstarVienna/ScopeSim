# -*- coding: utf-8 -*-
"""Contains tests for Quantization class."""

import pytest
import numpy as np

from scopesim.detector import Detector
from scopesim.effects.electronic import Quantization

from scopesim.tests.mocks.py_objects.header_objects import _implane_header


@pytest.fixture
def mock_detector():
    det = Detector(_implane_header())
    return det


@pytest.fixture(scope="function", name="det")
def detector_with_data(mock_detector):
    det = mock_detector
    width = det._hdu.data.shape[1]
    det._hdu.data[:] = 1.2
    det._hdu.data[:, width//2] = 1.99
    return det


class TestInit:
    def test_initialised_with_nothing(self):
        quant = Quantization()
        assert isinstance(quant, Quantization)


class TestApplyTo:
    def test_floors_pixels_to_integer_with_default(self, det):
        assert det.data.sum() > det.data.size * 1.2
        quant = Quantization()
        det = quant.apply_to(det)
        assert det.data.sum() == det.data.size
        assert det.data.dtype == np.uint32

    def test_floors_pixels_to_integer_with_custom(self, det):
        assert det.data.sum() > det.data.size * 1.2
        quant = Quantization(dtype="int16")
        det = quant.apply_to(det)
        assert det.data.sum() == det.data.size
        assert det.data.dtype == np.int16

    def test_logs_warning_for_float_dtype(self, det, caplog):
        assert det.data.sum() > det.data.size * 1.2
        quant = Quantization(dtype=float)
        det = quant.apply_to(det)
        assert det.data.sum() == det.data.size
        # Test sub-dtype because exact dtype for "float" may depend on OS
        assert np.issubdtype(det.data.dtype, np.floating)
        assert "Setting quantized data to dtype" in caplog.text
