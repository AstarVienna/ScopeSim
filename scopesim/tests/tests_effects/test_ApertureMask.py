import pytest
from pytest import approx

import numpy as np
from astropy.io import fits

from scopesim import rc
from scopesim.effects import ApertureMask
from scopesim.effects.apertures import points_on_a_circle, make_aperture_polygon
from scopesim.optics.fov_manager import FovVolumeList

import matplotlib.pyplot as plt


PLOTS = False


def basic_aperture_mask(x=(-2, -1, 1, 2), y=(-1, -2, 2, 1), **kwargs):
    params = {"array_dict": {"x": list(x), "y": list(y)},
              "x_unit": "arcsec",
              "y_unit": "arcsec"}
    params.update(kwargs)
    return ApertureMask(**params)


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(ApertureMask(), ApertureMask)

    def test_initialises_with_x_y_array_dict(self):
        apm = basic_aperture_mask()
        assert isinstance(apm, ApertureMask)
        assert "x" in apm.table.colnames

    def test_initialises_from_file(self, mock_path):
        apm = ApertureMask(filename=str(mock_path / "test_aperture.dat"))
        assert isinstance(apm, ApertureMask)
        assert "y" in apm.table.colnames


class TestApplyTo:
    def test_shrinks_fov_volume_list_as_expected(self):
        apm = basic_aperture_mask(x=[-2, 2, 2, -2], y=[-1, -1, 1, 1],
                                  pixel_scale=0.1)
        fvl = FovVolumeList()
        fvl = apm.apply_to(fvl)
        assert fvl.volumes[0]["x_min"] == -2

    def test_shrinks_fov_volume_list_with_extended_region(self):
        apm = basic_aperture_mask(x=[-2, 2, 2, -2], y=[-1, -1, 1, 1],
                                  extend_fov_beyond_slit=2, pixel_scale=0.1)
        fvl = FovVolumeList()
        fvl = apm.apply_to(fvl)
        assert fvl.volumes[0]["x_min"] == -4


class TestHeader:
    def test_returns_header_surrounding_tilted_aperture(self):
        e = 1e-7
        apm = basic_aperture_mask(x=[-2, -1, 1+e, 2+e],
                                  y=[-1, -2, 2+e, 1+e],
                                  pixel_scale=0.1)
        hdr = apm.header
        assert isinstance(hdr, fits.Header)
        assert hdr["NAXIS1"] == hdr["NAXIS2"] == 40


class TestMask:
    def test_returns_mask_that_covers_everything_for_zero_angle(self):
        kwargs = {"array_dict": {"x": [-1.5, 1.5, 1.5, -1.5],
                                 "y": [-0.01, -0.01, 0.01, 0.01]},
                  "x_unit": "arcsec",
                  "y_unit": "arcsec",
                  "pixel_scale": 0.004,
                  "no_mask": False}
        apm = ApertureMask(**kwargs)
        # known issue - for super thin apertures, the first row is masked
        # assert np.all(apm.mask)

        if PLOTS:
            plt.imshow(apm.mask.T)
            plt.show()

    def test_returns_mask_for_super_wierd_aperture(self):
        kwargs = {"array_dict": {"x": [-1, 1, 0, 1, -1, -0.5],
                                 "y": [-0.7, -1, 0, 1, 0.8, 0]},
                  "x_unit": "arcsec",
                  "y_unit": "arcsec",
                  "pixel_scale": 0.1,
                  "no_mask": False}
        apm = ApertureMask(**kwargs)

        if PLOTS:
            plt.imshow(apm.mask.T)
            plt.show()


class TestMakeAperturePolygon:
    @pytest.mark.parametrize("shape, n_corners", [("rect", 4), ("hex", 6),
                                                  ("round", 32), (7, 7)])
    def test_make_square(self, shape, n_corners):
        poly = make_aperture_polygon(-2, 2, 0.2, -0.2, 0, shape)
        assert len(poly["x"]) == n_corners
        assert np.max(poly["x"]) == approx(2, rel=1e-5)

        if PLOTS:
            plt.figure(figsize=(7, 7))
            plt.plot(poly["x"], poly["y"])
            plt.show()


class TestPointsOnACircle:
    def test_prints_a_circle(self):
        x, y = points_on_a_circle(32)
        assert np.max(x) == 1 and np.min(x) == -1

        if PLOTS:
            plt.figure(figsize=(7, 7))
            plt.plot(x, y)
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.show()
