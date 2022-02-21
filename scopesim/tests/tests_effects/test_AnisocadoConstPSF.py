import os
import pytest
from pytest import approx

import numpy as np
from astropy import units as u
from astropy.io import fits

from scopesim import rc
from scopesim.effects import FieldConstantPSF, AnisocadoConstPSF
from scopesim.tests.mocks.py_objects import fov_objects as fovobj
from scopesim.tests.mocks.py_objects import source_objects as srcobj

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


@pytest.fixture(scope="function")
def psf_object():
    psf = AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                            strehl=0.5, wavelength=2.15)
    return psf


@pytest.fixture(scope="function")
def fov_object():
    return fovobj._centre_micado_fov()


class TestInit:
    def test_throws_error_with_no_keywords(self):
        with pytest.raises(ValueError):
            AnisocadoConstPSF()

    def test_initialises_with_correct_input(self):
        psf = AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                                strehl=0.85, wavelength=2.15)
        assert isinstance(psf, AnisocadoConstPSF)

    def test_throws_error_if_desired_strehl_too_high(self):
        with pytest.raises(ValueError):
            AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                              strehl=0.99, wavelength=0.8)


@pytest.mark.usefixtures("psf_object", "fov_object")
class TestGetKernel:
    def test_strehl_map_is_in_data(self, psf_object):
        assert isinstance(psf_object._file[0], fits.PrimaryHDU)

    def test_returns_kernel(self, psf_object, fov_object):
        kernel = psf_object.get_kernel(fov_object)

        if PLOTS:
            plt.imshow(kernel, norm=LogNorm())
            plt.show()

        assert isinstance(kernel, np.ndarray)
        assert np.shape(kernel) == (512, 512)
        assert psf_object.strehl_ratio == approx(0.5, rel=0.01)

    def test_returns_kernel_for_filtername_wavelength(self):
        psf = AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                                strehl=0.15, wavelength="J")
        kernel = psf.get_kernel(0.004)

        if PLOTS:
            plt.imshow(kernel, norm=LogNorm())
            plt.show()

        assert isinstance(kernel, np.ndarray)


class TestApplyTo:
    def test_is_applied_to_point_sources(self):
        n = 10
        x, y, mag = 2 * np.random.random(size=(3, n)) - 1
        src = srcobj._vega_source(x=x[0], y=y[0], mag=mag[0])
        for i in range(1, n):
            src += srcobj._vega_source(x=x[i], y=y[i], mag=mag[i])
        fov = fovobj._centre_micado_fov(n=1)
        fov.extract_from(src)
        fov.view()

        psf = AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                                strehl=0.5, wavelength=2.15,
                                convolve_mode="same", psf_side_length=512)
        psf.apply_to(fov)

        if PLOTS:
            plt.imshow(fov.data, norm=LogNorm(), vmin=1)
            plt.show()

        assert 1e-99 < np.average(fov.data) < 1e99

