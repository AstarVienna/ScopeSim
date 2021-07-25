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
    kernel = np.zeros((n - 1, n - 1))
    kernel[:, n // 2] += 1
    kernel[n // 2, :] += 1
    kernel[n // 8 * 3:n // 8 * 5, n // 8 * 3:n // 8 * 5] += 1
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

