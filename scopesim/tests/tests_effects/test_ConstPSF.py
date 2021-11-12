# 1. ConstPSF should return the PSF for a given lambda
# 2. should throw errors when:
#   - file doesn't exist
#   - file doesn't have wave0 in the extension headers
# 3. should have attributes:
#   - lam_bin_centers : pulled from the header of the hduNs
#   - shape : (N_EXT, 1, NAXIS2, NAXIS1)
#   - nearest(wave, pos=None, hdu=False) : should return the array for the
#                                          given position and wavelength
#   - array : returns an array for the defaults values of wave and pos
#   - psf : returns self, as array returns a layer based on defaults
# 4. should have the following methods:
#   - (?) apply_to(obj)
#   - fov_grid(which="waveset")
#   - get_kernel(obj)

import os
import pytest
from pytest import approx

import numpy as np
from astropy import units as u
from astropy.io import fits

from scopesim import rc
from scopesim.optics.fov import FieldOfView
from scopesim.optics import image_plane_utils as imp_utils
from scopesim.effects import FieldConstantPSF, psfs

from scopesim.tests.mocks.py_objects.fov_objects import _centre_fov

import matplotlib.pyplot as plt

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


@pytest.fixture(scope="function")
def centre_fov():
    return _centre_fov()


class TestInit:
    def test_errors_when_initialised_with_nothing(self):
        with pytest.raises(ValueError):
            isinstance(FieldConstantPSF(), FieldConstantPSF)

    def test_initialised_when_passed_fits_filename(self):
        constpsf = FieldConstantPSF(filename="test_ConstPSF.fits")
        assert isinstance(constpsf, FieldConstantPSF)

    def test_initialised_when_passed_fits_filename_with_waveleng(self):
        constpsf = FieldConstantPSF(filename="test_ConstPSF_WAVELENG.fits",
                                    wave_key="WAVELENG")
        assert isinstance(constpsf, FieldConstantPSF)
        assert len(constpsf._waveset) == 2


class TestGetKernel:
    @pytest.mark.parametrize("waves, max_pixel",
                             [([1.1, 1.3], 1 / 3.),
                              ([1.5, 1.7], 1 / 3.),
                              ([1.9, 2.5], 1.),
                              ([0.9, 1.1], 1 / 3.),
                              ([1.875, 1.876], 1.)])
    def test_returns_array_with_single_kernel_from_fov(self, waves, max_pixel):
        constpsf = FieldConstantPSF(filename="test_ConstPSF.fits")
        fov = _centre_fov(n=10, waverange=waves)
        kernel = constpsf.get_kernel(fov)

        assert np.sum(kernel) == approx(1)
        assert np.max(kernel) == approx(max_pixel)

    @pytest.mark.parametrize("factor", [1/3., 1, 5])
    def test_kernel_is_scale_properly_if_cdelts_differ(self, factor):
        fov = _centre_fov(n=10, waverange=[1.5, 1.7])
        fov.header["CDELT1"] *= factor
        fov.header["CDELT2"] *= factor

        constpsf = FieldConstantPSF(filename="test_ConstPSF.fits")
        kernel = constpsf.get_kernel(fov)

        psf_shape = np.array(constpsf._file[2].data.shape)
        kernel_shape = kernel.shape

        assert np.all(kernel_shape == psf_shape / factor)


class TestApplyTo:
    @pytest.mark.parametrize("waves, max_pixel",
                             [([1.1, 1.3], 1 / 3.),
                              ([1.5, 1.7], 1 / 3.),
                              ([1.9, 2.5], 1.)])
    def test_convolves_with_basic_fov_for_each_waveleng(self, waves, max_pixel):
        centre_fov = _centre_fov(n=10, waverange=waves)
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.view("image")
        centre_fov.hdu.data = np.zeros((nax2, nax1))
        centre_fov.hdu.data[1::4, 1::4] = 1

        constpsf = FieldConstantPSF(filename="test_ConstPSF.fits")
        fov_returned = constpsf.apply_to(centre_fov)

        if PLOTS:
            plt.imshow(fov_returned.hdu.data, origin="lower")
            plt.show()

        assert np.max(fov_returned.hdu.data) == approx(max_pixel)

    def test_convolution_leaves_constant_background_intact(self):
        """FOV object with constant data must be unchanged by PSF convolution

        1. Set up initial FOV object with constant data array
        2. Convolve with PSF
        3. Assert that convolved array is identical to initial array.
        """
        centre_fov = _centre_fov(n=10, waverange=[1.1, 1.3])
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]

        centre_fov.view()
        centre_fov.hdu.data = np.ones((nax2, nax1), dtype=np.float32)

        constpsf = FieldConstantPSF(filename="test_ConstPSF.fits")
        fov_returned = constpsf.apply_to(centre_fov)

        if PLOTS:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7,3))
            plot_A = ax1.imshow(centre_fov.hdu.data)
            ax1.set_title("before convolution")
            fig.colorbar(plot_A, ax=ax1)
            plot_B = ax2.imshow(fov_returned.hdu.data)
            ax2.set_title("after convolution")
            fig.colorbar(plot_B, ax=ax2)
            plt.show()

        assert np.all(np.equal(fov_returned.hdu.data,
                               centre_fov.hdu.data))
