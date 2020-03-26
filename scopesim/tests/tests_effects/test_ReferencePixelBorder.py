import pytest
from pytest import approx

import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u

from scopesim.effects import electronic as ee
from scopesim.optics.image_plane import ImagePlane
from scopesim.base_classes import ImagePlaneBase

from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square

PLOTS = False


class TestInit:
    def test_initialised_with_nothing(self):
        rpb = ee.ReferencePixelBorder()
        assert isinstance(rpb, ee.ReferencePixelBorder)
        assert np.all([rpb.meta[key] == 0
                       for key in ["top", "bottom", "right", "left"]])

    def test_borders_all_set_to_5_for_keyword_all(self):
        rpb = ee.ReferencePixelBorder(all=5)
        assert np.all([rpb.meta[key] == 5
                       for key in ["top", "bottom", "right", "left"]])

    def test_border_set_differently(self):
        rpb = ee.ReferencePixelBorder(top=5, right=3)
        borders = {"top": 5, "bottom": 0, "right": 3, "left": 0}
        assert np.all([rpb.meta[key] == borders[key] for key in borders])


class TestApplyTo:
    def test_no_border_if_nothing_passed(self):
        implane = ImagePlane(_image_hdu_square().header)
        implane.hdu.data = np.ones(implane.hdu.data.shape)
        rpb = ee.ReferencePixelBorder()
        implane = rpb.apply_to(implane)

        assert np.sum(implane.data) == 10201

    def test_sets_border_to_zero(self):
        implane = ImagePlane(_image_hdu_square().header)
        implane.hdu.data = np.ones(implane.hdu.data.shape)
        rpb = ee.ReferencePixelBorder(all=5, top=15)
        implane = rpb.apply_to(implane)

        if PLOTS:
            plt.imshow(implane.data, origin="bottom")
            plt.show()

        assert np.sum(implane.data) == 7371
