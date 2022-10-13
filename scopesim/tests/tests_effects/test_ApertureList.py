import os
import pytest
from pytest import approx

import numpy as np
from astropy.io import fits

from scopesim import rc
from scopesim.effects import ApertureList, ApertureMask
from scopesim.optics.fov_manager import FovVolumeList

import matplotlib.pyplot as plt
PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
if FILES_PATH not in rc.__search_path__:
    rc.__search_path__ += [FILES_PATH]


def basic_aperture_list(**kwargs):
    params = {"array_dict": {"id": [0, 1, 2],
                             "left": [-3, -1, -1],
                             "right": [3, 1, 1],
                             "top": [-1, 1, 3],
                             "bottom": [-3, -1, 1],
                             "angle": [0, 0, 0],
                             "conserve_image": [True]*3,
                             "shape": ["rect", "hex", 7]},
              "x_unit": "arcsec",
              "y_unit": "arcsec",
              "pixel_scale": 0.01}
    params.update(kwargs)
    return ApertureList(**params)


class TestInit:
    def test_initialises_with_nothing(self):
        isinstance(ApertureList(), ApertureList)

    def test_initialises_with_array_dict(self):
        apl = basic_aperture_list()
        assert isinstance(apl, ApertureList)

    def test_initialises_with_filename(self):
        apl = ApertureList(filename="test_aperture_list.dat")
        assert isinstance(apl, ApertureList)


class TestApplyTo:
    def test_extracts_volumes_for_all_apertures(self):
        apl = basic_aperture_list()
        fvl = FovVolumeList()
        fvl = apl.apply_to(fvl)
        assert len(fvl.volumes) == 3
        assert fvl.volumes[0]["x_min"] == -3
        assert fvl.volumes[1]["x_min"] == -1

    def test_extracts_single_volume_for_all_apertures(self):
        apl = basic_aperture_list(fov_for_each_aperture=False)
        fvl = FovVolumeList()
        fvl = apl.apply_to(fvl)
        assert len(fvl.volumes) == 1
        assert fvl.volumes[0]["x_min"] == -3
        assert fvl.volumes[0]["y_max"] == 3

    def test_increases_all_volume_sizes_for_extend_beyond_slit(self):
        apl = basic_aperture_list(fov_for_each_aperture=True,
                                  extend_fov_beyond_slit=2)
        fvl = FovVolumeList()
        fvl = apl.apply_to(fvl)
        assert len(fvl.volumes) == 3
        assert fvl.volumes[0]["x_min"] == -5
        assert fvl.volumes[1]["x_min"] == -3

    def test_increases_single_volume_for_extend_beyond_slit(self):
        apl = basic_aperture_list(fov_for_each_aperture=False,
                                  extend_fov_beyond_slit=3)
        fvl = FovVolumeList()
        fvl = apl.apply_to(fvl)
        assert len(fvl.volumes) == 1
        assert fvl.volumes[0]["x_min"] == -6
        assert fvl.volumes[0]["y_max"] == 6


class TestApertures:
    def test_returns_list_of_aperture_masks(self):
        apl = ApertureList(filename="test_aperture_list.dat",
                           no_mask=False, pixel_scale=0.01)
        apertures = apl.apertures
        assert all([isinstance(am, ApertureMask) for am in apertures])

        if PLOTS:
            for ii in range(len(apertures)):
                plt.subplot(2, 2, ii+1)
                plt.imshow(apertures[ii].mask.T)
            plt.show()
