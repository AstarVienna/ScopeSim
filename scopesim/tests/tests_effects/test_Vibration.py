import pytest
from pytest import approx
import numpy as np
from astropy import units as u

from scopesim.effects import Vibration
from scopesim.optics.fov import FieldOfView
from scopesim.optics.image_plane import ImagePlane
from scopesim.tests.mocks.py_objects.header_objects import _fov_header, \
                                                           _implane_header

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


@pytest.fixture(scope="function")
def fov_hdr():
    return _fov_header()


@pytest.fixture(scope="function")
def implane_hdr():
    return _implane_header()


class TestInit:
    def test_initialises_with_required_parameters(self):
        kwargs = {"fwhm": 0.01, "pixel_scale": 0.004}
        assert isinstance(Vibration(**kwargs), Vibration)

    def test_throws_and_error_when_keyword_missing(self):
        kwargs = {"fwhm": 0.01}
        with pytest.raises(ValueError):
            Vibration(**kwargs)

    def test_kernel_sums_to_one(self):
        vibration = Vibration(**{"fwhm": 0.01, "pixel_scale": 0.004})
        vibration.get_kernel(None)
        assert np.sum(vibration.kernel) == approx(1, rel=1e-3)


@pytest.mark.usefixtures("fov_hdr", "implane_hdr")
class TestApplyTo:
    def test_nothing_happens_if_apply_to_fov(self, fov_hdr):
        fov = FieldOfView(header=fov_hdr, waverange=[0.5, 2.5], area=1*u.m**2)
        fov.view()
        fov.hdu.data = np.zeros((11, 11))
        fov.hdu.data[5, 5] = 1
        vibration = Vibration(**{"fwhm": 0.01, "pixel_scale": 0.004})
        vibration.apply_to(fov)

        assert fov.hdu.data[5, 5] == 1

    def test_something_happens_if_apply_to_imageplane(self, implane_hdr):
        implane = ImagePlane(header=implane_hdr)
        implane.hdu.data = np.zeros((11, 11))
        implane.hdu.data[5, 5] = 1
        vibration = Vibration(**{"fwhm": 0.2, "pixel_scale": 0.004})
        vibration.apply_to(implane)

        assert implane.hdu.data[5, 5] < 1

        if PLOTS:
            plt.imshow(vibration.kernel, norm=LogNorm())
            plt.colorbar()
            plt.show()


