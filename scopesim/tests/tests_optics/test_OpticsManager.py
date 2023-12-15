import pytest
from unittest.mock import patch
from astropy.io import fits

from scopesim.optics import optics_manager as opt_mgr
from scopesim.optics.optical_element import OpticalElement
from scopesim.effects import Effect

from scopesim.tests.mocks.py_objects.yaml_objects import\
    _inst_yaml_dict, _detector_yaml_dict


@pytest.fixture(scope="class")
def paths_patch(mock_path, mock_path_micado):
    with patch("scopesim.rc.__search_path__", [mock_path, mock_path_micado]):
        yield


@pytest.fixture(scope="function")
def inst_yaml_dict():
    return _inst_yaml_dict()


@pytest.fixture(scope="function")
def detector_yaml_dict():
    return _detector_yaml_dict()


@pytest.mark.usefixtures("paths_patch")
class TestOpticsManager:
    def test_initialises_with_nothing(self):
        assert isinstance(opt_mgr.OpticsManager(),
                          opt_mgr.OpticsManager)

    def test_initialises_yaml_dict(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager(detector_yaml_dict)
        assert isinstance(opt_man, opt_mgr.OpticsManager)

    def test_initialises_two_yaml_dicts(self, detector_yaml_dict,
                                        inst_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict, inst_yaml_dict])
        assert isinstance(opt_man, opt_mgr.OpticsManager)
        assert len(opt_man.optical_elements) == 2

    def test_has_effects_loaded(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        # print(opt_man.optical_elements[1])
        assert isinstance(opt_man.optical_elements[0],
                          opt_mgr.OpticalElement)
        assert isinstance(opt_man.optical_elements[0].effects[0], Effect)


@pytest.mark.usefixtures("patch_mock_path")
class TestOpticsManagerImagePlaneHeader:
    def test_makes_image_plane_header_correctly(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager(detector_yaml_dict)
        opt_man.meta["SIM_PIXEL_SCALE"] = 0.004
        print(opt_man)
        assert isinstance(opt_man.image_plane_headers[0], fits.Header)


@pytest.mark.usefixtures("patch_mock_path")
class TestGetItem:
    def test_returns_optical_element(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        assert isinstance(opt_man["micado_detector_array"], OpticalElement)

    def test_returns_effect(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        assert isinstance(opt_man["detector_qe_curve"], Effect)

    def test_raise_error_for_wrong_name(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        with pytest.raises(ValueError):
            opt_man["knut_dietrich"]

    def test_returns_value_from_simple_hash_string(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        value = opt_man["#detector_qe_curve.filename"]
        assert value == "TER_blank.dat"

    def test_returns_value_from_full_hash_string(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        value = opt_man["#micado_detector_array.detector_qe_curve.filename"]
        assert value == "TER_blank.dat"

    @pytest.mark.parametrize("key", [("detector_qe_curve.filename"),
                                     ("#detector_qe_curve")])
    def test_errors_on_wrong_hash_strings(self, detector_yaml_dict, key):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        with pytest.raises(ValueError):
            opt_man[key]
