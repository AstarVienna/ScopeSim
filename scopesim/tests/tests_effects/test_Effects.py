import os

import numpy as np
import pytest
from astropy.table import Table

import scopesim as sim
from scopesim.effects import Effect
from scopesim.effects import ApertureList
from scopesim.effects import SurfaceList

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__ += [MOCK_PATH]


@pytest.fixture()
def surf_list_file():
    fname = os.path.join(MOCK_PATH, "LIST_mirrors_MICADO_Wide.tbl")
    return fname


@pytest.fixture()
def aperture_list_file():
    fname = os.path.join(MOCK_PATH, "LIST_mirrors_MICADO_Wide.tbl")
    return fname


class TestEffectInit:
    def test_initialises_with_no_input(self):
        assert isinstance(Effect(), Effect)

    def test_initalising_with_arrays_creates_table(self):
        eff = Effect(array_dict={"x": [-1, 0, 1], "y": [1, 0, 1],
                                 "flux": [1, 2, 3]})
        assert isinstance(eff, Effect)
        assert np.sum(eff.table["flux"]) == 6

    def test_has_method_apply_to(self):
        assert hasattr(Effect(), "apply_to")

    def test_has_method_fov_grid(self):
        assert hasattr(Effect(), "fov_grid")


@pytest.mark.usefixtures("surf_list_file")
class TestSurfaceListInit:
    def test_initialises_with_nothing(self):
        assert isinstance(SurfaceList(), SurfaceList)

    def test_initialises_with_valid_filename(self, surf_list_file):
        surf_list = SurfaceList(filename=surf_list_file)
        assert isinstance(surf_list, SurfaceList)
        assert isinstance(surf_list.data, Table)


@pytest.mark.usefixtures("aperture_list_file")
class TestApertureListInit:
    def test_initialises_with_nothing(self):
        assert isinstance(ApertureList(), ApertureList)

    def test_initialises_with_valid_filename(self, surf_list_file):
        aperture_list = ApertureList(filename=surf_list_file)
        assert isinstance(aperture_list, ApertureList)
        assert isinstance(aperture_list.data, Table)
