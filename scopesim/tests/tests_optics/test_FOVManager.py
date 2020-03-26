import os
import pytest
from pytest import approx

import numpy as np

import scopesim as sim
from scopesim import rc
from scopesim.optics.fov_manager import FOVManager
from scopesim.optics import fov_manager as fov_mgr
from scopesim.optics.image_plane import ImagePlane
from scopesim.commands import UserCommands

from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim.tests.mocks.py_objects.yaml_objects import \
     _usr_cmds_min_viable_scope

import matplotlib.pyplot as plt

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


@pytest.fixture(scope="function")
def mvs_effects_list():
    return eo._mvs_effects_list()


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
                x, y = calc_footprint(fov.hdu.header)
                plt.fill(x*3600, y*3600)

            plt.subplot(122)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.hdu.header, "D")
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


@pytest.mark.usefixtures("ap_list", "spt_list", "det_list",
                         "config_yaml", "spec_source")
class TestGenerateFOVsSpectroscopyMode:
    def test_input_data_checks_out(self, ap_list, spt_list, det_list,
                                   config_yaml):
        if PLOTS:
            plt.subplot(121)
            for hdr in hdrs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(hdr)
                plt.plot(x*3600, y*3600)
                plt.title("Sky plane")
                plt.xlabel("[arcsec]")

            plt.subplot(122)
            for hdr in hdrs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(hdr, "D")
                plt.plot(x, y)
                plt.title("Detector focal plane")
                plt.xlabel("[mm]")

            plt.show()


class TestGetImagingFOVs:
    def test_returns_FOV_objects_for_basic_input(self):
        apm = eo._img_aperture_mask(array_dict={"x": [-1.0, 1.0, 1.0, -1.0],
                                                "y": [-1.0, -1.0, 1.0, 1.0]})
        kwargs = {"pixel_scale": 0.01, "plate_scale": 1,
                  "max_segment_size": 100 ** 2, "chunk_size": 100}

        hdrs = fov_mgr.get_imaging_headers([apm], **kwargs)
        waveset = np.linspace(1, 2, 6)
        shifts = {"wavelengths": np.array([1, 2]),
                  "x_shifts": np.zeros(2),
                  "y_shifts": np.array([0, 1]) / 3600}  # 0..1 arcsec shift
        fovs = fov_mgr.get_imaging_fovs(headers=hdrs, waveset=waveset,
                                        shifts=shifts)

        assert len(fovs) == (len(waveset)-1) * len(hdrs)

        if PLOTS:
            from scopesim.optics.image_plane_utils import calc_footprint
            plt.subplot(121)
            for fov in fovs:
                x, y = calc_footprint(fov.hdu.header)
                plt.fill(x*3600, y*3600, alpha=0.1, c="b")
                plt.title("Sky plane")
                plt.xlabel("[arcsec]")

        if PLOTS:
            det_arr = DetectorArray(det_list)
            hdus = det_arr.readout([implane])

            for i, hdu in enumerate(hdus[1:5]):
                plt.subplot(2, 2, i+1)
                plt.imshow(hdu.data, origin="lower")
            plt.show()


class TestGetSpectroscopyFovs2:
    # fov_mgr.get_spectroscopy_fovs2()
    pass

