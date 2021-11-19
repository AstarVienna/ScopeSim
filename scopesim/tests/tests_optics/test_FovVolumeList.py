import pytest
from pytest import approx
import numpy as np
from astropy import units as u
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.optics.fov_manager import FovVolumeList


class TestInit:
    def test_init_with_nothing(self):
        fvl = FovVolumeList()
        assert isinstance(fvl, FovVolumeList)
        assert len(fvl.volumes) == 1

class TestSplit:
    @pytest.mark.parametrize("axis, value", [("wave", 3.0),
                                             ("x", 0.),
                                             ("y", -900)])
    def test_split_inside_limits_makes_two_volumes(self, axis, value):
        fvl = FovVolumeList()
        fvl.split(axis=axis, value=value)

        assert len(fvl) == 2
        assert fvl[0][f"{axis}_max"] == value
        assert fvl[1][f"{axis}_min"] == value

    def test_split_twice_along_single_axis_makes_three_volumes(self):
        fvl = FovVolumeList()
        fvl.split(axis="wave", value=[3.0, 3.1])

        assert len(fvl) == 3
        for i, value in enumerate([3.0, 3.1, 30.]):
            assert fvl[i]["wave_max"] == value

    def test_split_along_three_dimensions_makes_eight_volumes(self):
        fvl = FovVolumeList()
        fvl.split(axis=["wave", "x", "y"], value=[3.0, 0, -900])

        assert len(fvl) == 8

    def test_two_splits_along_three_dimensions_makes_27_volumes(self):
        fvl = FovVolumeList()
        fvl.split(axis=["wave", "x", "y"], value=([3.0, 3.1], [-1, 1], [0, 5]))

        assert len(fvl) == 27


class TestShrink:
    def test_shrinks_limits_for_single_volume(self):
        fvl = FovVolumeList()
        fvl.shrink("wave", [2.9, 3.1])

        assert len(fvl) == 1
        assert fvl[0]["wave_min"] == 2.9

    def test_shrinks_limits_only_for_relevent_volume(self):
        fvl = FovVolumeList()
        fvl.split(axis="wave", value=[3.0, 3.1])
        fvl.shrink("wave", [2.9, 3.2])

        assert len(fvl) == 3
        assert fvl[0]["wave_min"] == 2.9
        assert fvl[1]["wave_min"] == 3.0
        assert fvl[2]["wave_min"] == 3.1
        assert fvl[2]["wave_max"] == 3.2

    def test_doesnt_shrink_if_values_too_large(self):
        fvl = FovVolumeList()
        fvl.shrink("wave", [0.1, 35])

        assert fvl[0]["wave_min"] == 0.3
        assert fvl[0]["wave_max"] == 30

    def test_shrink_along_two_axes(self):
        fvl = FovVolumeList()
        fvl.shrink(["x", "y"], ([0.1, None], [None, 5]))
        print(fvl)

        assert fvl[0]["x_min"] == 0.1
        assert fvl[0]["y_max"] == 5
        assert fvl[0]["y_min"] == -1800

    def test_removes_volumes_no_longer_inside_volume_limits(self):
        fvl = FovVolumeList()
        fvl.split("wave", [3.0, 3.1])

        assert len(fvl) == 3

        fvl.shrink("wave", [3.0, 3.1000001])

        assert len(fvl) == 2


class TestExtract:
    def test_extracts_volume_fully_inside_old_volume(self):
        fvl = FovVolumeList()
        new_vols = fvl.extract(["x", "y"], ([-1, 1], [0, 5]))

        assert new_vols[0]["x_min"] == -1
        assert new_vols[0]["y_max"] == 5

    def test_extracts_volume_partially_inside_old_volume(self):
        fvl = FovVolumeList()
        fvl.shrink("x", [0, None])
        new_vols = fvl.extract(["x", "y"], ([-1, 1], [0, 5]))

        assert new_vols[0]["x_min"] == 0
        assert new_vols[0]["x_max"] == 1

    def test_extracts_no_volume_outside_old_volume(self):
        fvl = FovVolumeList()
        fvl.shrink("x", [2, None])
        new_vols = fvl.extract(["x", "y"], ([-1, 1], [0, 5]))

        assert len(new_vols) == 0

    def test_extracts_two_volumes_inside_split_old_volume(self):
        fvl = FovVolumeList()
        fvl.split("x", 0)
        new_vols = fvl.extract(["x", "y"], ([-1, 1], [0, 5]))

        assert len(new_vols) == 2

    def test_extracts_only_wavelength_range(self):
        fvl = FovVolumeList()
        fvl.split("x", 0)
        new_vols = fvl.extract(["wave"], ([1, 2], ))

        assert len(new_vols) == 2


class TestAdd:
    def test_adds_extracted_volumes(self):
        fvl = FovVolumeList()
        fvl.split("x", 0)
        new_vols = fvl.extract(["x", "y"], ([-1, 1], [0, 5]))
        fvl = fvl + new_vols

        assert len(fvl) == 4

    def test_iadds_extracted_volumes(self):
        fvl = FovVolumeList()
        fvl.split("x", 0)
        new_vols = fvl.extract(["x", "y"], ([-1, 1], [0, 5]))
        fvl += new_vols

        assert len(fvl) == 4
