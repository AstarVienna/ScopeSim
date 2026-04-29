"""
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
"""

import pytest
from pytest import approx

import numpy as np
from numpy import testing as npt

import matplotlib.pyplot as plt

from scopesim.effects import FieldConstantPSF

from scopesim.tests.mocks.py_objects.fov_objects import _centre_fov



PLOTS = False

# pylint: disable=missing-class-docstring,
# pylint: disable=missing-function-docstring

@pytest.fixture(scope="function")
def centre_fov():
    return _centre_fov()


@pytest.fixture(scope="function")
def const_psf(mock_path):
    constpsf = FieldConstantPSF(
        filename=str(mock_path / "test_ConstPSF.fits"))
    return constpsf


class TestInit:
    def test_errors_when_initialised_with_nothing(self):
        with pytest.raises(ValueError):
            isinstance(FieldConstantPSF(), FieldConstantPSF)

    def test_initialised_when_passed_fits_filename(self, const_psf):
        assert isinstance(const_psf, FieldConstantPSF)

    def test_initialised_when_passed_fits_filename_with_waveleng(self,
                                                                 mock_path):
        constpsf = FieldConstantPSF(
            filename=str(mock_path / "test_ConstPSF_WAVELENG.fits"),
            wave_key="WAVELENG")
        assert isinstance(constpsf, FieldConstantPSF)
        assert len(constpsf._waveset) == 2

    def test_initialised_with_psf_name_and_file_format(self, mock_path):
        # actually tests code in parent class DiscretePSF
        psf_name = "ConstPSF"
        file_format = "test_{}_WAVELENG.fits"
        constpsf = FieldConstantPSF(
            psf_name=psf_name,
            psf_path=mock_path,
            filename_format=file_format,
        )
        assert isinstance(constpsf, FieldConstantPSF)


class TestGetKernel:
    @pytest.mark.parametrize("waves, max_pixel",
                             [([1.1, 1.3], 1 / 3.),
                              ([1.5, 1.7], 1 / 3.),
                              ([1.9, 2.5], 1.),
                              ([0.9, 1.1], 1 / 3.),
                              ([1.875, 1.876], 1.)])
    def test_returns_array_with_single_kernel_from_fov(
        self,
        waves,
        max_pixel,
        const_psf,
    ):
        fov = _centre_fov(n=10, waverange=waves)
        fov.view()
        kernel = const_psf.get_kernel(fov)

        npt.assert_allclose(kernel.sum(), 1)
        npt.assert_allclose(kernel.max(), max_pixel)

    @pytest.mark.parametrize("scale_factor", (
        pytest.param(1/3, marks=pytest.mark.xfail(reason="PSF center off")),
        1,
        pytest.param(5, marks=pytest.mark.xfail(reason="1x1 array turned into scalar")),
    ))
    def test_kernel_is_scale_properly_if_cdelts_differ(
        self,
        scale_factor,
        const_psf,
    ):
        fov = _centre_fov(n=10, waverange=[1.5, 1.7])
        fov.header["CDELT1"] *= scale_factor
        fov.header["CDELT2"] *= scale_factor
        fov.view()

        kernel = const_psf.get_kernel(fov)
        kernel_shape = np.array(kernel.shape)

        # Check if kernel is scaled by correct factor in both dimensions
        npt.assert_allclose(
            const_psf._file[2].data.shape / kernel_shape,
            scale_factor,
        )
        # Check if max of kernel is in the center pixel
        assert (tuple(kernel_shape // 2) ==
                np.unravel_index(kernel.argmax(), kernel.shape))


class TestApplyTo:
    @pytest.mark.parametrize("waves, max_pixel",
                             [([1.1, 1.3], 1 / 3.),
                              ([1.5, 1.7], 1 / 3.),
                              ([1.9, 2.5], 1.)])
    def test_convolves_with_basic_fov_for_each_waveleng(
        self,
        waves,
        max_pixel,
        const_psf,
    ):
        centre_fov = _centre_fov(n=10, waverange=waves)
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.view("image")
        centre_fov.hdu.data = np.zeros((nax2, nax1))
        centre_fov.hdu.data[1::4, 1::4] = 1

        fov_returned = const_psf.apply_to(centre_fov)

        if PLOTS:
            plt.imshow(fov_returned.hdu.data, origin="lower")
            plt.show()

        npt.assert_allclose(fov_returned.hdu.data.max(), max_pixel)

    def test_convolution_leaves_constant_background_intact(self, const_psf):
        """FOV object with constant data must be unchanged by PSF convolution

        1. Set up initial FOV object with constant data array
        2. Convolve with PSF
        3. Assert that convolved array is identical to initial array.
        """
        centre_fov = _centre_fov(n=10, waverange=[1.1, 1.3])
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]

        centre_fov.view()
        centre_fov.hdu.data = np.ones((nax2, nax1), dtype=np.float32)

        fov_returned = const_psf.apply_to(centre_fov)

        if PLOTS:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7,3))
            plot_a = ax1.imshow(centre_fov.hdu.data)
            ax1.set_title("before convolution")
            fig.colorbar(plot_a, ax=ax1)
            plot_a = ax2.imshow(fov_returned.hdu.data)
            ax2.set_title("after convolution")
            fig.colorbar(plot_a, ax=ax2)
            plt.show()

        npt.assert_array_equal(fov_returned.hdu.data, centre_fov.hdu.data)
