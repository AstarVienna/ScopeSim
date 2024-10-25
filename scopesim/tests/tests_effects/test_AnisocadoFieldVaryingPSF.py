import pytest
from pytest import approx

import numpy as np
from astropy.io import fits

from scopesim.effects.psfs import semianalytical as sa
from scopesim.tests.mocks.py_objects import fov_objects as fovobj
from scopesim.tests.mocks.py_objects import source_objects as srcobj

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


PLOTS = False


def mock_psf_cube(n=64):
    psf_cube = sa.generate_anisocado_psf_cube(dxs=[0, 0, 10, 10],
                                              dys=[0, 10, 0, 10],  # arcsec
                                              n=n,
                                              pixel_scale=0.01,
                                              wave=3.4)
    return psf_cube


class TestMakeXys:
    def test_returns_coords_in_square_mode(self, r_max=3, n=4):
        xs, ys = sa.make_xys(r_max, n, grid_type="square")

        assert len(xs) == n**2
        assert xs[-1] == ys[-1] == r_max
        assert xs[0] == ys[0] == 0

    def test_returns_coords_in_radial_mode(self, r_max=3, n=4):
        xs, ys = sa.make_xys(r_max, n, grid_type="radial")

        assert len(xs) == len(ys) == n * (n-1) + 1
        assert xs[-1] == approx(0)
        assert ys[-1] == r_max
        assert xs[0] == ys[0] == 0


class TestMakeDecimalIndexMap:
    def test_returns_imagehdu_with_useful_header_info(self, n=64):
        xs, ys = sa.make_xys(r_max=9, n=4, grid_type="radial")
        idx_im = sa.make_decimal_index_map(xs, ys, n)

        if PLOTS:
            plt.imshow(idx_im.data, origin="lower")
            plt.show()

        assert isinstance(idx_im, fits.ImageHDU)
        assert idx_im.data.shape == (2*n, 2*n)
        assert idx_im.data[0, 1] - 0.75 == idx_im.data[-1, -2]


class TestGenerateAnisocadoPsfCube:
    def test_makes_cube_of_expected_shape(self, n=32):
        psf_cube = mock_psf_cube(n).data
        assert psf_cube.shape == (4, n, n)


class TestGetPsfFromDecimalIndex:
    def test_returns_image_for_sensible_input(self, n=32):
        psf_cube = mock_psf_cube(n).data
        psf_im = sa.get_psf_from_decimal_index(psf_cube, 1.25)

        if PLOTS:
            plt.imshow(psf_im, norm="log")
            plt.show()

        assert psf_im.shape == (n, n)
