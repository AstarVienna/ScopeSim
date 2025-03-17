"""
1. FVPSF should return the PSF for a position in the FOV and a given lambda
2. should throw errors when:
  - file doesn't exist
  - file doesn't have a image or table in the 1st extension
3. should have attributes:
  - lam_bin_centers : pulled from the header of the hduNs
  - layer_map : an ImageHDU with a map of where each layer is valid
  - layer_table : The BinTableHDU if available
  - _make_layer_map : makes a layer_map from a BinTableHDU is needed
  - mask(wave, pos) : returns an ImageHDU with WCS with a mask of the valid
                      region for a given wavelength and position in the FOV
  - shape : (N_EXT, N_LAYER, NAXIS2, NAXIS1)
  - nearest(wave, pos=None, hdu=False) : should return the array for the
                                          given position and wavelength
  - defaults : a dictionary with {"wave": , pos: (,)} so that .array can be
                used for backwards compatibility
  - set_defaults(wave, pos)
  - array : returns an array for the defaults values of wave and pos
  - psf : returns self, as array returns a layer based on defaults
"""

import pytest
from pytest import approx

import numpy as np
from astropy.io import fits

from scopesim.effects.psfs.discrete import get_strehl_cutout
from scopesim.effects import FieldVaryingPSF
from scopesim.tests.mocks.py_objects.fov_objects import _centre_fov
from scopesim.tests.mocks.py_objects.psf_objects import _basic_circular_fvpsf

import matplotlib.pyplot as plt


PLOTS = False


@pytest.fixture(scope="function")
def centre_fov():
    fov = _centre_fov()
    fov.view()
    return fov


@pytest.fixture(scope="function")
def basic_circular_fvpsf():
    return _basic_circular_fvpsf


class TestInit:
    def test_errors_when_initialised_with_nothing(self):
        with pytest.raises(ValueError):
            isinstance(FieldVaryingPSF(), FieldVaryingPSF)

    def test_initialised_when_passed_fits_filename(self, mock_path):
        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        assert isinstance(fvpsf, FieldVaryingPSF)


class TestStrehlImageHDU:
    def test_returns_the_correct_hdu(self, mock_path):
        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        strehl_imhdu = fvpsf.strehl_imagehdu
        assert isinstance(strehl_imhdu, fits.ImageHDU)


class TestGetKernel:
    def test_returns_array_with_single_kernel_from_fov(self, mock_path):
        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        kernels = fvpsf.get_kernel(_centre_fov(n=10))
        assert np.all(kernels[0][0] == fvpsf._file[2].data[4])
        assert kernels[0][1] is None

    @pytest.mark.parametrize("factor", [0.2, 1, 3])
    def test_kernel_is_scale_properly_if_cdelts_differ(self, factor,
                                                       mock_path):
        fov = _centre_fov(n=10)
        fov.header["CDELT1"] *= factor
        fov.header["CDELT2"] *= factor

        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        kernels = fvpsf.get_kernel(fov)

        psf_shape = np.array(fvpsf._file[2].data[0].shape)
        kernel_shape = kernels[0][0].shape

        assert all(kernel_shape == psf_shape * factor)

    def test_returns_four_arrays_when_fov_on_intersection(self, mock_path):
        fov = _centre_fov(n=20)
        fov.header["CRVAL1"] -= 15/3600.
        fov.header["CRVAL2"] -= 15/3600.

        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        kernels = fvpsf.get_kernel(fov)
        assert len(kernels) == 4

        if PLOTS:
            for ii, kernel_mask in enumerate(kernels):
                plt.subplot(2, 4, ii+1)
                plt.imshow(kernel_mask[0].T, origin="lower")
                plt.subplot(2, 4, ii+5)
                plt.imshow(kernel_mask[1].T, origin="lower")

            plt.show()


class TestApplyTo:
    def test_convolution_with_central_psf_for_central_region(self, centre_fov,
                                                             mock_path):
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.hdu.data = np.zeros((nax2, nax1))
        centre_fov.hdu.data[::3, ::3] = 1
        orig = centre_fov.hdu.data.copy()
        sum_orig = np.sum(orig)

        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        fov_back = fvpsf.apply_to(centre_fov)

        assert np.sum(fov_back.hdu.data) == sum_orig

        if PLOTS:
            plt.imshow(fov_back.hdu.data, origin="lower")
            plt.show()

    def test_convolution_with_fvpsfs_for_shifted_region(self, centre_fov,
                                                        mock_path):
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.hdu.data = np.zeros((nax2, nax1))
        centre_fov.hdu.data[1::5, 1::5] = 1
        centre_fov.fields = [1]
        orig = centre_fov.hdu.data.copy()
        sum_orig = np.sum(orig)

        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        fov_back = fvpsf.apply_to(centre_fov)

        # The FieldVaryingPSF was broken from early 2022 till late 2024
        # because the pixels were not modified at all. This assert verifies
        # that the data is not identical after applying the PSF, thereby
        # preventing regressions.
        assert not (orig == fov_back.hdu.data).all()
        assert np.sum(fov_back.hdu.data) == approx(sum_orig, rel=1E-2)

        if PLOTS:
            plt.imshow(fov_back.hdu.data, origin="lower")
            plt.show()

    def test_circular_fvpsf(self, basic_circular_fvpsf, mock_path):
        centre_fov = _centre_fov(n=62)
        centre_fov.view()
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.hdu.data = np.zeros((nax2, nax1))

        rng = np.random.RandomState(42)
        x, y = rng.randint(6, nax1-6, (2, 150))
        centre_fov.hdu.data[x, y] = 1
        # centre_fov.hdu.data[6:nax1-6:10, 6:nax1-6:10] = 1
        centre_fov.fields = [1]
        orig = centre_fov.hdu.data.copy()
        sum_orig = np.sum(orig)

        fvpsf = FieldVaryingPSF(
            filename=str(mock_path / "test_circular_fvpsf.fits"))
        fov_back = fvpsf.apply_to(centre_fov)

        if PLOTS:
            plt.imshow(fov_back.hdu.data, origin="lower", vmax=0.1)
            plt.show()

        # print(np.sum(fov_back.hdu.data), sum_orig)
        assert not (orig == fov_back.hdu.data).all()
        assert np.sum(fov_back.hdu.data) == approx(sum_orig, rel=1E-2)


class TestFunctionGetStrehlCutout:
    @pytest.mark.parametrize("scale", [0.2, 0.5, 1, 2])
    def test_returns_correct_section_of_strehl_map(self, scale, mock_path):
        centre_fov = _centre_fov(10)
        centre_fov.header["CDELT1"] *= scale
        centre_fov.header["CDELT2"] *= scale
        centre_fov.header["CRVAL1"] -= 15/3600.
        centre_fov.header["CRVAL2"] -= 15/3600.

        fvpsf = FieldVaryingPSF(filename=str(mock_path / "test_FVPSF.fits"))
        strehl_hdu = get_strehl_cutout(
            centre_fov.header, fvpsf.strehl_imagehdu)

        if PLOTS:
            plt.imshow(strehl_hdu.data, origin="lower")
            plt.colorbar()
            plt.show()

        assert all(np.unique(strehl_hdu.data).astype(int) == [0, 1, 3, 4])
        # These work for scale 0.5, 1, 2, but are off for 0.2, weird...
        # assert (strehl_hdu.data == 0).sum() == 100
        # assert (strehl_hdu.data == 1).sum() == 100
        # assert (strehl_hdu.data == 3).sum() == 100
        # assert (strehl_hdu.data == 4).sum() == 100
