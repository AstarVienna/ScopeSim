import pytest

from scopesim.optics import optical_element as opt_elem
from scopesim.optics.effects import GaussianDiffractionPSF
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
