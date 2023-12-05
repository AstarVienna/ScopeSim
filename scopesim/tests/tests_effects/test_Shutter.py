# -*- coding: utf-8 -*-
"""Contains tests for Shutter class."""

import pytest
import numpy as np

from scopesim.effects.shutter import Shutter
from scopesim.optics.image_plane import ImagePlane

from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square


@pytest.fixture(scope="function", name="implane")
def mock_implane():
    implane = ImagePlane(_image_hdu_square().header)
    implane.hdu.data = np.ones(implane.hdu.data.shape)
    return implane


class TestInit:
    def test_initialised_with_nothing(self):
        shtr = Shutter()
        assert isinstance(shtr, Shutter)


class TestApplyTo:
    def test_sets_pixels_to_zero(self, implane):
        assert implane.data.sum() == implane.data.size
        shtr = Shutter()
        implane = shtr.apply_to(implane)
        assert implane.data.sum() == 0.0

    def test_shape_is_preserved(self, implane):
        orig_shape = implane.data.shape
        shtr = Shutter()
        implane = shtr.apply_to(implane)
        assert implane.data.shape == orig_shape
