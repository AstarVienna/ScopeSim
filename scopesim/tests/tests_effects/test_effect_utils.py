import pytest

from scopesim import effects as efs
from scopesim.effects import effects_utils as eu, GaussianDiffractionPSF, \
    SurfaceList
from scopesim.tests.mocks.py_objects.effects_objects import _surf_list, \
    _surf_list_empty, _filter_surface
from scopesim.tests.mocks.py_objects.yaml_objects import _atmo_yaml_dict


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


class TestMakeEffect:
    def test_it_creates_an_effects_object(self, atmo_yaml_dict):
        effdic = atmo_yaml_dict["effects"][0]
        properties = atmo_yaml_dict["properties"]
        effect = eu.make_effect(effdic, **properties)

        assert isinstance(effect, GaussianDiffractionPSF)
        assert effect.meta["diameter"] == 39


class TestCombineSurfaceEffects:
    def test_load_just_one_surface(self, filter_surface):
        surf_list = eu.combine_surface_effects([filter_surface])
        assert isinstance(surf_list, SurfaceList)
        assert len(surf_list.table) == 1

    def test_load_two_surfaces(self, filter_surface):
        surf_list = eu.combine_surface_effects([filter_surface] * 3)
        assert len(surf_list.table) == 3

    def test_load_just_one_surface_list(self, surf_list):
        new_surf_list = eu.combine_surface_effects([surf_list])
        assert len(new_surf_list.table) == len(surf_list.table)

    def test_load_two_surface_lists(self, surf_list):
        new_surf_list = eu.combine_surface_effects([surf_list] * 3)
        assert len(new_surf_list.table) == 3 * len(surf_list.table)

    def test_load_one_surface_and_one_surface_list(self, surf_list,
                                                   filter_surface):
        new_surf_list = eu.combine_surface_effects([filter_surface, surf_list])
        assert len(new_surf_list.table) == len(surf_list.table) + 1


class TestScopesimEffectClasses:
    def test_all_effects_including_effect_base(self):
        all_efs = eu.scopesim_effect_classes()
        assert efs.Effect in list(all_efs.values())
        assert len(all_efs) > 2

    def test_only_psf_effects_returned(self):
        all_efs = eu.scopesim_effect_classes(efs.PSF)
        assert efs.Effect not in list(all_efs.values())
        assert all(["psf" in eff for eff in all_efs])
