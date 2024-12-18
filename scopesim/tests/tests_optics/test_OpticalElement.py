import pytest
from unittest.mock import patch

from scopesim.optics import optical_element as opt_elem
from scopesim.effects import GaussianDiffractionPSF
from scopesim.commands import UserCommands
from scopesim.effects import Effect

from scopesim.tests.mocks.py_objects.yaml_objects import _atmo_yaml_dict, \
    _detector_yaml_dict


@pytest.fixture(scope="function")
def atmo_yaml_dict():
    return _atmo_yaml_dict()


@pytest.fixture(scope="function")
def detector_yaml_dict():
    return _detector_yaml_dict()


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

    # TODO: Check whether the test below actually does what I think it's testing.
    # Should be obsolete once the currsys is removed from rc
    def test_currsys_ignore_effects_have_false_include_flag(self, atmo_yaml_dict):
        with patch("scopesim.rc.__currsys__", UserCommands()) as patched:
            patched.ignore_effects = ["super_psf"]
            opt_el = opt_elem.OpticalElement(atmo_yaml_dict, cmds=patched)
            for ii in patched.ignore_effects:
                assert opt_el["super_psf"].include is False


@pytest.mark.usefixtures("patch_mock_path")
class TestOpticalElementGetZOrderEffects:
    @pytest.mark.parametrize("z_lvl, zmax, n", [(0, None, 2),
                                                (100, None, 1),
                                                (200, 299, 1)])
    def test_returns_the_effects_with_z_values(self, z_lvl, zmax, n,
                                               detector_yaml_dict):
        opt_el = opt_elem.OpticalElement(detector_yaml_dict)
        assert len(list(opt_el.get_z_order_effects(z_lvl, zmax))) == n


@pytest.mark.usefixtures("patch_mock_path")
class TestGetItem:
    def test_returns_effect_for_normal_string(self, detector_yaml_dict):
        opt_el = opt_elem.OpticalElement(detector_yaml_dict)
        eff = opt_el["detector_qe_curve"]
        assert isinstance(eff, Effect)

    def test_returns_effect_for_hash_string(self, detector_yaml_dict):
        opt_el = opt_elem.OpticalElement(detector_yaml_dict)
        value = opt_el["#detector_qe_curve.filename"]
        assert value == "TER_blank.dat"

    @pytest.mark.parametrize("key", [("detector_qe_curve.filename"),
                                     ("#detector_qe_curve")])
    def test_pass_silently_when_not_valid_hash_string(self, detector_yaml_dict,
                                                     key):
        opt_el = opt_elem.OpticalElement(detector_yaml_dict)
        assert len(opt_el[key]) == 0


class TestOpticalElementSurfaceListProperty:
    def test_returns_empty_list_if_no_surface_list_given(self):
        pass
