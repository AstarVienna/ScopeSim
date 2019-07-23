import os
import scopesim as sim

import pytest
from astropy.io import fits

from scopesim.optics import optics_manager as opt_mgr
from scopesim.effects import Effect

from scopesim.tests.mocks.py_objects.yaml_objects import\
    _inst_yaml_dict, _detector_yaml_dict


FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__ += [FILES_PATH, YAMLS_PATH]


@pytest.fixture(scope="function")
def inst_yaml_dict():
    return _inst_yaml_dict()


@pytest.fixture(scope="function")
def detector_yaml_dict():
    return _detector_yaml_dict()


@pytest.mark.usefixtures("detector_yaml_dict", "inst_yaml_dict")
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
        assert len(opt_man.optical_elements) == 3

    def test_has_effects_loaded(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager([detector_yaml_dict])
        # print(opt_man.optical_elements[1])
        assert isinstance(opt_man.optical_elements[1],
                          opt_mgr.OpticalElement)
        assert isinstance(opt_man.optical_elements[1].effects[0], Effect)


@pytest.mark.usefixtures("detector_yaml_dict")
class TestOpticsManagerImagePlaneHeader:
    def test_makes_image_plane_header_correctly(self, detector_yaml_dict):
        opt_man = opt_mgr.OpticsManager(detector_yaml_dict)
        opt_man.meta["SIM_PIXEL_SCALE"] = 0.004
        print(opt_man)
        assert isinstance(opt_man.image_plane_header, fits.Header)
