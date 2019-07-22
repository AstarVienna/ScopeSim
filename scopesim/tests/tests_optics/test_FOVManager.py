import os
import pytest

import numpy as np

import scopesim as sim
from scopesim.optics.fov_manager import FOVManager
from scopesim.optics import fov_manager as fov_mgr
from scopesim.optics.image_plane import ImagePlane

from scopesim.tests.mocks.py_objects.effects_objects import _mvs_effects_list
from scopesim.tests.mocks.py_objects.yaml_objects import \
     _usr_cmds_min_viable_scope

import matplotlib.pyplot as plt

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
sim.rc.__search_path__ += [FILES_PATH, YAMLS_PATH]


@pytest.fixture(scope="function")
def mvs_effects_list():
    return _mvs_effects_list()


@pytest.fixture(scope="function")
def mvs_usr_cmds():
    return _usr_cmds_min_viable_scope()


@pytest.mark.usefixtures("mvs_effects_list")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(), FOVManager)

    def test_initialises_with_list_of_effects(self, mvs_effects_list):
        assert isinstance(FOVManager(mvs_effects_list), FOVManager)


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGenerateFovs:
    def test_returns_the_desired_number_of_fovs(self, mvs_effects_list,
                                                mvs_usr_cmds):
        for yaml_dic in mvs_usr_cmds:
            sim.rc.__currsys__.update(yaml_dic)
        fov_man = FOVManager(mvs_effects_list)
        fovs = fov_man.generate_fovs_list()

        implane_size = mvs_effects_list[-2].image_plane_header["NAXIS1"]
        chunk_size = sim.rc.__currsys__["!SIM.computing.chunk_size"]
        wave_layers = 1

        assert len(fovs) == np.round(implane_size / chunk_size)**2 * wave_layers

        if PLOTS:
            plt.subplot(121)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.hdu.header)
                plt.plot(x, y)

            plt.subplot(122)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.hdu.header, "D")
                plt.plot(x, y)

            plt.show()

    def test_fovs_dont_overlap_on_canvas(self, mvs_effects_list, mvs_usr_cmds):

        for yaml_dic in mvs_usr_cmds:
            sim.rc.__currsys__.update(yaml_dic)

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


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGet3DShifts:
    def test_returns_empty_for_the_moment(self):
        pass

    def test_returns_zeros_when_adc_and_ad_are_equal(self):
        pass

    def test_returns_offset_when_adc_has_set_airmass(self):
        pass


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGetImagingWaveset:
    def test_returns_waveset_from_cmds(self):
        pass

    def test_returns_waveset_based_on_psfs(self):
        pass


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGetImagingHeaders:
    def test_returns_header_grid_detector_size(self):
        pass


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGetImagingFOVs:
    def test_returns_header_grid_detector_size(self):
        pass

