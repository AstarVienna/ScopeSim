import os
import pytest

from scopesim import rc
from scopesim.effects import effects_utils as e_utils, GaussianDiffractionPSF, \
    SurfaceList
from scopesim.tests.mocks.py_objects.effects_objects import _surf_list, \
    _surf_list_empty, _filter_surface
from scopesim.tests.mocks.py_objects.yaml_objects import _atmo_yaml_dict

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]


@pytest.fixture(scope="function")
def atmo_yaml_dict():
    return _atmo_yaml_dict()


@pytest.fixture(scope="function")
def surf_list():
    return _surf_list()


@pytest.fixture(scope="function")
def surf_list_empty():
    return _surf_list_empty()


@pytest.fixture(scope="function")
def filter_surface():
    return _filter_surface()


@pytest.mark.usefixtures("atmo_yaml_dict")
class TestMakeEffect:
    def test_it_creates_an_effects_object(self, atmo_yaml_dict):
        effdic = atmo_yaml_dict["effects"][0]
        properties = atmo_yaml_dict["properties"]
        effect = e_utils.make_effect(effdic, **properties)

        assert isinstance(effect, GaussianDiffractionPSF)
        assert effect.meta["diameter"] == 39


@pytest.mark.usefixtures("surf_list", "filter_surface")
class TestMakeRadiometryTable:
    def test_load_just_one_surface(self, filter_surface):
        rad_table = e_utils.combine_surface_effects([filter_surface])
        assert isinstance(rad_table, SurfaceList)

    def test_load_two_surfaces(self, filter_surface):
        rad_table = e_utils.combine_surface_effects([filter_surface] * 3)
        assert len(rad_table.radiometry_table.table) == 3

    def test_load_just_one_surface_list(self, surf_list):
        rad_table = e_utils.combine_surface_effects([surf_list])
        len1 = len(surf_list.table)
        len2 = len(rad_table.radiometry_table.table)
        assert len2 == len1

    def test_load_two_surface_lists(self, surf_list):
        rad_table = e_utils.combine_surface_effects([surf_list, surf_list])
        len1 = len(surf_list.radiometry_table.table)
        len2 = len(rad_table.radiometry_table.table)
        assert len2 == 2 * len1

    def test_load_one_surface_and_one_surface_list(self, surf_list,
                                                   filter_surface):
        rad_table = e_utils.combine_surface_effects([surf_list,
                                                     filter_surface])
        len1 = len(surf_list.table)
        len2 = len(rad_table.table)
        assert len2 == len1 + 1

