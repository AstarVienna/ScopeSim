import pytest

from scopesim import rc
from scopesim.optics import optical_element as opt_elem
from scopesim.effects import GaussianDiffractionPSF
from scopesim.commands import UserCommands

from scopesim.tests.mocks.py_objects.yaml_objects import _atmo_yaml_dict, \
    _detector_yaml_dict


@pytest.fixture(scope="function")
def atmo_yaml_dict():
    return _atmo_yaml_dict()


@pytest.fixture(scope="function")
def detector_yaml_dict():
    return _detector_yaml_dict()


@pytest.mark.usefixtures("atmo_yaml_dict")
class TestOpticalElementInit:
    def test_initialised_with_nothing(self):
        assert isinstance(opt_elem.OpticalElement(),
                          opt_elem.OpticalElement)

    def test_initialised_with_yaml_dict(self, atmo_yaml_dict):
        opt_el = opt_elem.OpticalElement(atmo_yaml_dict)
        assert isinstance(opt_el, opt_elem.OpticalElement)
        assert isinstance(opt_el.effects[0], GaussianDiffractionPSF)

    def test_cleans_obs_keywords_from_yaml_properties(self, atmo_yaml_dict):
        atmo_yaml_dict["properties"]["tempertaure"] = "OBS_TEMPERATURE"
        obs_dict = {"OBS_TEMPERATURE": 0}
        opt_el = opt_elem.OpticalElement(atmo_yaml_dict, **obs_dict)

        assert opt_el.properties["temperature"] == 0

    def test_make_not_included_effects_but_set_flag_false(self, atmo_yaml_dict):
        opt_el = opt_elem.OpticalElement(atmo_yaml_dict)
        assert len(opt_el.effects) == 3
        assert opt_el.effects[2].include is False

    def test_currsys_ignore_effects_have_false_include_flag(self, atmo_yaml_dict):
        rc.__currsys__ = UserCommands()
        rc.__currsys__.ignore_effects = ["super_psf", "atmo_dispersion"]
        opt_el = opt_elem.OpticalElement(atmo_yaml_dict)
        for ii in range(2):
            assert opt_el.effects[ii].include is False


@pytest.mark.usefixtures("detector_yaml_dict")
class TestOpticalElementGetZOrderEffects:
    @pytest.mark.parametrize("z_orders, n", [(0, 2), (100, 1), ([200, 299], 1)])
    def test_returns_the_effects_with_z_values(self, z_orders, n,
                                               detector_yaml_dict):
        opt_el = opt_elem.OpticalElement(detector_yaml_dict)
        assert len(opt_el.get_z_order_effects(z_orders)) == n


@pytest.mark.usefixtures("detector_yaml_dict")
class TestOpticalElementSurfaceListProperty:
    def test_returns_empty_list_if_no_surface_list_given(self):
        pass

