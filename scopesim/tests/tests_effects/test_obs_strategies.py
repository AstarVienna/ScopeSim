import pytest

import numpy as np
import matplotlib.pyplot as plt

from scopesim.effects.obs_strategies import ChopNodCombiner
from scopesim.detector import Detector
from scopesim.tests.mocks.py_objects.header_objects import _implane_header


PLOTS = False


@pytest.fixture(scope="function")
def basic_detector():
    # An image plane image with background=100 and object_flux=1
    det = Detector(_implane_header())
    n = det.header["NAXIS1"]
    im = np.ones((n, n)) * 100              # background level = 100
    im[n//2-5:n//2+5, n//2-5:n//2+5] += 1   # object flux level += 1
    det.hdu.data = im

    return det


class TestChopNodCombinerInit:
    def test_initilises_with_required_keywords(self):
        cnc = ChopNodCombiner(pixel_scale=0.004, chop_offsets=(0.12, 0))
        assert isinstance(cnc, ChopNodCombiner)

    def test_throws_error_when_missing_keywords(self):
        with pytest.raises(ValueError):
            ChopNodCombiner()


class TestChopNodCombinerApplyTo:
    def test_creates_image_for_parallel_chop(self, basic_detector):
        cnc = ChopNodCombiner(pixel_scale=0.004, chop_offsets=(0.12, 0))
        basic_detector = cnc.apply_to(basic_detector)

        if PLOTS:
            plt.imshow(basic_detector.hdu.data)
            plt.show()

        outimg = basic_detector.hdu.data
        assert outimg.sum() == 0
        assert ((outimg == 0.).sum() / outimg.size) > .8  # most elements zero

    def test_creates_image_for_perpendicular_chop(self, basic_detector):
        cnc = ChopNodCombiner(pixel_scale=0.004, chop_offsets=(0.12, 0),
                              nod_offsets=(0, -0.20))
        basic_detector = cnc.apply_to(basic_detector)

        if PLOTS:
            plt.imshow(basic_detector.hdu.data)
            plt.show()

        outimg = basic_detector.hdu.data
        assert outimg.sum() == 0
        assert ((outimg == 0.).sum() / outimg.size) > .8  # most elements zero
