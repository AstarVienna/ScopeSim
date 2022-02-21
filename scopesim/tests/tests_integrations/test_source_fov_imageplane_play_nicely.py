# pylint: disable=no-self-use
# pylint: disable=missing-function-docstring
# pylint: disable=invalid-name

import numpy as np
import pytest
from pytest import approx

from astropy import units as u

from scopesim.optics.fov import FieldOfView
from scopesim.optics.image_plane import ImagePlane

from scopesim.tests.mocks.py_objects import source_objects as src
from scopesim.tests.mocks.py_objects import header_objects as hdrs

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


@pytest.fixture(scope="function")
def comb_src():
    return src._combined_source()


@pytest.fixture(scope="function")
def fov_hdr():
    return hdrs._fov_header()


@pytest.fixture(scope="function")
def implane_hdr():
    return hdrs._implane_header()


@pytest.mark.usefixtures("comb_src", "fov_hdr", "implane_hdr")
class TestInteractionBetweenSourceFOVImagePlane:
    """
    Test:
    - fov extracts correctly from source object and converts to image
    - fov image is correctly added to the image plane
    """
    def test_can_extract_the_source_in_a_fov(self, fov_hdr, comb_src,
                                             implane_hdr):

        fov = FieldOfView(fov_hdr, waverange=[0.5, 2.5]*u.um, area=1*u.m**2)
        imp = ImagePlane(implane_hdr)

        fov.extract_from(comb_src)
        fov.view()
        imp.add(fov.hdu, wcs_suffix="D")

        assert np.sum(imp.image) > 0

        if PLOTS:
            plt.subplot(131)
            plt.imshow(comb_src.fields[3].data, origin="lower", norm=LogNorm())

            plt.subplot(132)
            plt.imshow(fov.hdu.data, origin="lower", norm=LogNorm())

            plt.subplot(133)
            plt.imshow(imp.hdu.data, origin="lower", norm=LogNorm())
            plt.show()
