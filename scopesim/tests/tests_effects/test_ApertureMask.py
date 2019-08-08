import os
import pytest
from pytest import approx

import numpy as np
from astropy.io import fits

from scopesim import rc
from scopesim.effects import ApertureMask

import matplotlib.pyplot as plt
PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
if FILES_PATH not in rc.__search_path__:
    rc.__search_path__ += [FILES_PATH]


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(ApertureMask(), ApertureMask)

    def test_initialises_with_x_y_array_dict(self):
        kwargs = {"array_dict": {"x": [-2, -1, 1, 2],
                                 "y": [-1, -2, 2, 1]},
                  "x_unit": "arcsec",
                  "y_unit": "arcsec"}
        apm = ApertureMask(**kwargs)
        assert isinstance(apm, ApertureMask)
        assert "x" in apm.table.colnames

    def test_initialises_from_file(self):
        apm = ApertureMask(filename="test_aperture.dat")
        assert isinstance(apm, ApertureMask)
        assert "y" in apm.table.colnames


class TestHeader:
    def test_returns_header_surrounding_tilted_aperture(self):
        e = 1e-7
        kwargs = {"array_dict": {"x": [-2, -1, 1+e, 2+e],
                                 "y": [-1, -2, 2+e, 1+e]},
                  "x_unit": "arcsec",
                  "y_unit": "arcsec",
                  "pixel_scale": 0.1}
        apm = ApertureMask(**kwargs)
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
        assert np.all(apm.mask)

        if PLOTS:
            plt.imshow(apm.mask.T)
            plt.show()


class TestFovGrid:
    # not needed because all it does is return the header
    pass
