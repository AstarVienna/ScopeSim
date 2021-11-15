import pytest
from pytest import approx
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim.effects import PSF
from scopesim.effects.psf_utils import rotational_blur
from scopesim.optics import ImagePlane
from scopesim.tests.mocks.py_objects.header_objects import _implane_header


PLOTS = False


def basic_image_plane():
    implane = ImagePlane(_implane_header())

    return implane


def basic_kernel(n=128):
    kernel = np.zeros((n, n))
    kernel[:, n // 2] += 1
    kernel[n // 2, :] += 1
    w0, w1 = int(np.floor(n / 8 * 3)), int(np.ceil(n / 8 * 5))
    kernel[w0:w1, w0:w1] += 1
    kernel /= np.sum(kernel)

    return kernel


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(PSF(), PSF)


class TestGetKernel:
    def test_returns_array_from_kernel(self):
        psf = PSF()
        psf.kernel = basic_kernel()

        assert np.sum(psf.kernel) == approx(np.sum(psf.get_kernel(None)))


class TestRotationBlur:
    def test_returns_rotated_kernel_array_has_same_sum(self):
        # Without blur
        implane = basic_image_plane()
        implane.data[75,75] = 1

        psf = PSF(rotational_blur_angle=0)
        psf.kernel = basic_kernel()
        implane = psf.apply_to(implane)

        if PLOTS:
            plt.subplot(121)
            plt.imshow(implane.data)

        # Without blur
        implane = basic_image_plane()
        implane.data[75,75] = 1

        psf = PSF(rotational_blur_angle=15)
        psf.kernel = basic_kernel()
        implane = psf.apply_to(implane)

        if PLOTS:
            plt.subplot(122)
            plt.imshow(implane.data)
            plt.show()

        assert np.sum(implane.data) == approx(1)


class TestApplyTo:
    def test_convolves_with_2D_image(self):
        implane = basic_image_plane()
        implane.data[75,75] = 1

        if not PLOTS:
            plt.subplot(131)
            plt.imshow(implane.data)

        psf = PSF()
        psf.kernel = basic_kernel(n=15)
        implane = psf.apply_to(implane)

        if PLOTS:
            plt.subplot(132)
            plt.imshow(psf.kernel)

            plt.subplot(133)
            plt.imshow(implane.data)
            plt.show()

        assert np.sum(implane.data) == approx(1)
        assert implane.data[75, 75] == approx(np.max(psf.kernel))

    def test_convolves_with_3D_cube(self):
        implane = basic_image_plane()
        implane.data[75,75] = 1
        implane.data = implane.data[None, :, :] * np.ones(3)[:, None, None]

        if not PLOTS:
            plt.subplot(231)
            plt.imshow(implane.data[1, :, :])

        psf = PSF(rotational_blur_angle=15, bkg_width=5)
        psf.kernel = basic_kernel(n=63)
        implane = psf.apply_to(implane)

        if PLOTS:
            plt.subplot(232)
            plt.imshow(psf.kernel)

            for i in range(3):
                plt.subplot(2,3,4+i)
                plt.imshow(implane.data[i, :, :])
            plt.show()

        assert np.sum(implane.data[1, :, :]) == approx(1)
        assert implane.data[1, 75, 75] == approx(np.max(psf.kernel), rel=1e-2)
