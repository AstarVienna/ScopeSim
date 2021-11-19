import os
import pytest

import numpy as np

import scopesim as sim
from scopesim import rc
from scopesim.optics import FieldOfView, FOVManager, ImagePlane
from scopesim.commands import UserCommands

from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim.tests.mocks.py_objects import yaml_objects as yo
from scopesim.tests.mocks.py_objects import integr_spectroscopy_objects as iso

import matplotlib.pyplot as plt

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


################################################################################
# Everything needed to test the FOVManager in Spectroscopy mode


@pytest.fixture(scope="function")
def ap_list():
    return iso.mock_aperture_list()


@pytest.fixture(scope="function")
def spt_list():
    return iso.mock_spectral_trace_list()


@pytest.fixture(scope="function")
def det_list():
    return iso.mock_detector_list()


@pytest.fixture(scope="function")
def config_yaml():
    return iso.mock_config_yaml()


@pytest.fixture(scope="function")
def spec_source():
    return iso.mock_point_source_object()


################################################################################
# Needed for FOVManager in Imaging mode

@pytest.fixture(scope="function")
def mvs_effects_list():
    return eo._mvs_effects_list()


@pytest.fixture(scope="function")
def mvs_usr_cmds():
    return yo._usr_cmds_min_viable_scope()


################################################################################


@pytest.mark.usefixtures("mvs_effects_list")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(preload_fovs=False), FOVManager)

    def test_initialises_with_list_of_effects(self, mvs_effects_list):
        assert isinstance(FOVManager(mvs_effects_list,
                                     preload_fovs=False), FOVManager)


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGenerateFovsImagingMode:
    def test_returns_the_desired_number_of_fovs(self, mvs_effects_list,
                                                mvs_usr_cmds):
        rc.__currsys__ = UserCommands(yamls=mvs_usr_cmds)
        rc.__currsys__["!SIM.computing.max_segment_size"] = 2**20

        fov_man = FOVManager(mvs_effects_list)
        fovs = fov_man.generate_fovs_list()

        implane_size = mvs_effects_list[-2].image_plane_header["NAXIS1"]
        chunk_size = sim.rc.__currsys__["!SIM.computing.chunk_size"]
        wave_layers = 1

        assert len(fovs) == np.round(implane_size / chunk_size)**2 * wave_layers

        # ..todo:: go searching for the white line
        if PLOTS:
            plt.subplot(121)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.header)
                plt.fill(x*3600, y*3600)

            plt.subplot(122)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.header, "D")
                plt.plot(x, y)

            plt.show()

    def test_fovs_dont_overlap_on_canvas(self, mvs_effects_list, mvs_usr_cmds):
        rc.__currsys__ = UserCommands(yamls=mvs_usr_cmds)

        implane = ImagePlane(mvs_effects_list[-2].image_plane_header)
        fov_man = FOVManager(mvs_effects_list)
        fovs = fov_man.generate_fovs_list()
        for fov in fovs:
            fov.hdu.data = np.ones((fov.header["NAXIS1"], fov.header["NAXIS2"]))
            implane.add(fov.hdu, wcs_suffix="D")

        if PLOTS:
            plt.imshow(implane.image.T, origin="lower")
            plt.colorbar()
            plt.show()

        # .. todo: Fix this - find out why!!!
        # assert np.all(implane.image == 1)

@pytest.mark.skip(reason="Ignoring old Spectroscopy integration tests")
@pytest.mark.usefixtures("ap_list", "spt_list", "det_list",
                         "config_yaml", "spec_source")
class TestGenerateFOVsSpectroscopyMode:
    def test_input_data_checks_out(self, ap_list, spt_list, det_list,
                                   config_yaml):
        if PLOTS:
            plt.subplot(121)
            ap_list.plot()
            plt.subplot(122)
            spt_list.plot(0.8, 2.5)
            det_list.plot()
            plt.show()

    def test_generates_fovs_for_spectroscopy_mode(self, ap_list, spt_list,
                                                  det_list, config_yaml,
                                                  spec_source):
        # Uses a basic SpectralTraceList with 4 traces and 2 apertures
        # 2 traces attached to each apertures
        config_yaml["wave_min"] = 1.0
        config_yaml["wave_max"] = 2.5

        effects = [ap_list, spt_list, det_list]
        fov_mgr = FOVManager(effects=effects, **config_yaml)
        fovs = fov_mgr.fovs
        print(len(fovs))

        assert all([isinstance(fov, FieldOfView) for fov in fovs])

        if PLOTS:
            implane = ImagePlane(det_list.image_plane_header)
            for fov in fovs:
                fov.extract_from(spec_source)
                fov.view()
                implane.add(fov.hdu, wcs_suffix="D")

            plt.imshow(implane.data, origin="lower")
            plt.show()
