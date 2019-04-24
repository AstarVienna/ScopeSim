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
sim.rc.__search_path__ += [FILES_PATH]


@pytest.fixture(scope="function")
def mvs_effects_list():
    return _mvs_effects_list()


@pytest.fixture(scope="function")
def mvs_usr_cmds():
    return _usr_cmds_min_viable_scope()


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(), FOVManager)

    def test_initialises_with_list_of_effects(self, mvs_effects_list):
        assert isinstance(FOVManager(mvs_effects_list), FOVManager)


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGenerateFovs:
    def test_returns_the_desired_number_of_fovs(self, mvs_effects_list,
                                                mvs_usr_cmds):
        fov_man = FOVManager(mvs_effects_list, **mvs_usr_cmds)
        fovs = fov_man.generate_fovs_list()

        implane_size = mvs_effects_list[-2].image_plane_header["NAXIS1"]
        chunk_size = mvs_usr_cmds["SIM_CHUNK_SIZE"]
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

        implane = ImagePlane(mvs_effects_list[-2].image_plane_header)
        fov_man = FOVManager(mvs_effects_list, **mvs_usr_cmds)
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
    def test_returns_empty_for_the_moment(self, mvs_effects_list, mvs_usr_cmds):
        shift_dict = fov_mgr.get_3d_shifts(mvs_effects_list, **mvs_usr_cmds)
        assert shift_dict["wavelengths"][0] == mvs_usr_cmds["SIM_LAM_MIN"]

    def test_returns_zeros_when_adc_and_ad_are_equal(self):
        pass

    def test_returns_offset_when_adc_has_set_airmass(self):
        pass


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGetImagingWaveset:
    def test_returns_waveset_from_cmds(self, mvs_effects_list, mvs_usr_cmds):
        waveset = fov_mgr.get_imaging_waveset(mvs_effects_list, **mvs_usr_cmds)
        assert np.all(waveset == [mvs_usr_cmds["SIM_LAM_MIN"],
                                  mvs_usr_cmds["SIM_LAM_MAX"]])

    def test_returns_waveset_based_on_psfs(self):
        pass


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGetImagingHeaders:
    def test_returns_header_grid_detector_size(self, mvs_effects_list,
                                               mvs_usr_cmds):
        hdrs = fov_mgr.get_imaging_headers(mvs_effects_list, **mvs_usr_cmds)
        implane_size = mvs_effects_list[-2].image_plane_header["NAXIS1"]
        chunk_size = mvs_usr_cmds["SIM_CHUNK_SIZE"]
        print(np.round(implane_size / chunk_size))

        assert len(hdrs) == np.round(implane_size / chunk_size)**2


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGetImagingFOVs:
    def test_returns_header_grid_detector_size(self, mvs_effects_list,
                                               mvs_usr_cmds):
        shift_dict = fov_mgr.get_3d_shifts(mvs_effects_list, **mvs_usr_cmds)
        waveset = fov_mgr.get_imaging_waveset(mvs_effects_list, **mvs_usr_cmds)
        hdrs = fov_mgr.get_imaging_headers(mvs_effects_list, **mvs_usr_cmds)

        hdrs = fov_mgr.get_imaging_fovs(hdrs, waveset, shift_dict)

        implane_size = mvs_effects_list[-2].image_plane_header["NAXIS1"]
        chunk_size = mvs_usr_cmds["SIM_CHUNK_SIZE"]
        wave_layers = len(waveset) - 1

        assert len(hdrs) == np.round(implane_size / chunk_size)**2 * wave_layers


