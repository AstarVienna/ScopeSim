import pytest
from pytest import approx

import numpy as np
from astropy import units as u

from scopesim import effects as efs
from scopesim.optics.fov import FieldOfView
from scopesim.optics.image_plane_utils import pix2val

from scopesim.tests.mocks.py_objects.source_objects import _image_source
from scopesim.tests.mocks.py_objects.header_objects import _basic_fov_header

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


def _basic_fov():
    src = _image_source()
    fov = FieldOfView(_basic_fov_header(), waverange=[1, 2]*u.um, area=1*u.m**2)
    fov.extract_from(src)

    return fov


@pytest.fixture(scope="function")
def basic_fov():
    return _basic_fov()


class TestInit:
    def test_does_not_initialise_with_nothing(self):
        with pytest.raises(TypeError):
            efs.GaussianDiffractionPSF()

    def test_initialises_with_only_diameter(self):
        eff = efs.GaussianDiffractionPSF(1)
        assert isinstance(eff, efs.GaussianDiffractionPSF)

    def test_initialised_with_other_keywords(self):
        eff = efs.GaussianDiffractionPSF(1, sub_pixel=False)
        assert eff.meta["sub_pixel"] is False


@pytest.mark.usefixtures("basic_fov")
class TestApplyTo:
    def test_size_of_fov_increases_when_convolved_if_convolve_mode_set_to_full(self, basic_fov):
        effect = efs.GaussianDiffractionPSF(1, convolve_mode="full")
        basic_fov.view()
        #basic_fov = effect.apply_to(basic_fov)

        orig_size = np.prod(basic_fov.fields[0].data.shape)
        orig_sum = np.sum(basic_fov.fields[0].data) * 2     # integ. flux of spec = 2.0
        new_size = np.prod(basic_fov.hdu.data.shape)
        new_sum = np.sum(basic_fov.hdu.data)

        assert new_size > orig_size
        assert new_sum == approx(orig_sum, rel=1e-3)

        if PLOTS:
            plt.subplot(121)
            plt.imshow(basic_fov.fields[0].data, origin="lower", norm=LogNorm())
            plt.subplot(122)
            plt.imshow(basic_fov.data, origin="lower", norm=LogNorm())
            plt.show()


    def test_crval_stays_the_same(self, basic_fov):
        basic_fov.header["CRPIX1"] = 0
        basic_fov.view()
        effect = efs.GaussianDiffractionPSF(1)

        x0, y0 = pix2val(basic_fov.header,
                         basic_fov.header["CRPIX1"],
                         basic_fov.header["CRPIX2"])

        basic_fov = effect.apply_to(basic_fov)
        x1, y1 = pix2val(basic_fov.header,
                         basic_fov.header["CRPIX1"],
                         basic_fov.header["CRPIX2"])

        assert x0 == x1 and y0 == y1
