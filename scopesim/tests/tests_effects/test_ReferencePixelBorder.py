"""Tests for ReferencePixelBorder"""
# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
import pytest
import numpy as np

from scopesim.effects import ReferencePixelBorder

from scopesim.tests.mocks.py_objects.detector_objects import _basic_detector


class TestInit:
    def test_initialised_with_nothing(self):
        rpb = ReferencePixelBorder()
        assert isinstance(rpb, ReferencePixelBorder)
        assert rpb.meta["border"] == [0, 0, 0, 0]

    def test_initialised_with_border(self):
        rpb = ReferencePixelBorder(border=[4, 5, 6, 7])
        assert rpb.meta["border"] == [4, 5, 6, 7]

    def test_error_with_short_border(self):
        with pytest.raises(ValueError):
            _ = ReferencePixelBorder(border=[1,3])

    def test_error_with_long_border(self):
        with pytest.raises(ValueError):
            _ = ReferencePixelBorder(border=[2, 3, 4, 5, 6])


class TestApplyTo:
    def test_no_border_if_nothing_passed(self):
        det = _basic_detector(width=128)
        det.hdu.data = np.ones(det.data.shape)
        rpb = ReferencePixelBorder()
        rpb.apply_to(det)
        assert np.sum(det.data) == 128**2

    def test_sets_border_to_zero(self):
        det = _basic_detector(width=128)
        det.hdu.data = np.ones(det.data.shape)
        rpb = ReferencePixelBorder(border=[32, 32, 32, 32])
        rpb.apply_to(det)
        assert np.sum(det.data) == 64**2

    def test_order_of_parameters_is_correct(self):
        det = _basic_detector(width=128)
        det.hdu.data = np.ones(det.data.shape)
        rpb = ReferencePixelBorder(border=[10, 20, 30, 0])
        rpb.apply_to(det)
        assert np.all(det.data[:10, ] == 0)
        assert np.all(det.data[10:(128-30), 20:] == 1)
        assert np.all(det.data[128-30+1:, :] == 0)
        assert np.all(det.data[:, :20] == 0)

    def test_all_zero_when_border_exeeds_image_limits(self):
        det = _basic_detector(width=128)
        det.hdu.data = np.ones(det.data.shape, dtype=np.float32)
        rpb = ReferencePixelBorder(border=[0, 400, 1000, 0])
        rpb.apply_to(det)
        assert np.all(det.data == 0)
