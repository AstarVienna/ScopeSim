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


class TestScopesimEffectClasses:
    def test_all_effects_including_effect_base(self):
        all_efs = eu.scopesim_effect_classes()
        assert efs.Effect in list(all_efs.values())
        assert len(all_efs) > 2

    def test_only_psf_effects_returned(self):
        all_efs = eu.scopesim_effect_classes(efs.PSF)
        assert efs.Effect not in list(all_efs.values())
        assert all(["psf" in eff for eff in all_efs])
