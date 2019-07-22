import os
from copy import deepcopy
import pytest
from pytest import approx

import numpy as np
from astropy import units as u

import scopesim as sim
from scopesim.optics.fov_manager import FOVManager
from scopesim.optics.image_plane import ImagePlane
from scopesim.optics.optical_train import OpticalTrain
from scopesim.optics.optics_manager import OpticsManager
from scopesim.commands.user_commands import UserCommands
from scopesim.utils import find_file

from scopesim.tests.mocks.py_objects import source_objects as src_objs

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))

sim.rc.__search_path__ += [FILES_PATH, YAMLS_PATH]


def _basic_cmds():
    return UserCommands(yamls=[find_file("CMD_mvs_cmds.yaml")])


def _unity_cmds():
    return UserCommands(yamls=[find_file("CMD_unity_cmds.yaml")])


@pytest.fixture(scope="function")
def cmds():
    return _basic_cmds()\


@pytest.fixture(scope="function")
def unity_cmds():
    return _unity_cmds()


@pytest.fixture(scope="function")
def tbl_src():
    return src_objs._table_source()


@pytest.fixture(scope="function")
def im_src():
    return src_objs._image_source()


@pytest.fixture(scope="function")
def unity_src():
    return src_objs._unity_source()


@pytest.mark.usefixtures("cmds")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(OpticalTrain(), OpticalTrain)

    def test_initialises_with_basic_commands(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt, OpticalTrain)

    def test_has_user_commands_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.cmds, UserCommands)

    def test_has_optics_manager_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.optics_manager, OpticsManager)

    def test_has_fov_manager_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        print(opt.fov_manager)
        assert isinstance(opt.fov_manager, FOVManager)

    def test_has_image_plane_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.image_plane, ImagePlane)

    def test_has_yaml_dict_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.yaml_dicts, list) and len(opt.yaml_dicts) > 0


@pytest.mark.usefixtures("cmds", "im_src", "tbl_src")
class TestObserve:
    """
    All tests here are for visual inspection.
    No asserts, this just to test that the puzzle gets put back together
    after it is chopped up by the FOVs.
    """

    def test_observe_works_for_table(self, cmds, tbl_src):
        opt = OpticalTrain(cmds)
        opt.observe(tbl_src)

        if PLOTS:
            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm())
            plt.show()

    def test_observe_works_for_image(self, cmds, im_src):
        opt = OpticalTrain(cmds)
        opt.observe(im_src)

        if PLOTS:
            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm())
            plt.show()

    def test_observe_works_for_source_distributed_over_several_fovs(self, cmds,
                                                                    im_src):

        orig_sum = np.sum(im_src.fields[0].data)

        cmds["SIM_PIXEL_SCALE"] = 0.02
        opt = OpticalTrain(cmds)
        opt.observe(im_src)

        wave = np.arange(0.5, 2.51, 0.1)*u.um
        unit = u.Unit("ph s-1 m-2 um-1")
        print(opt.optics_manager.surfaces_table.emission(wave).to(unit))
        print(opt.optics_manager.surfaces_table.table)
        final_sum = np.sum(opt.image_plane.image)
        print(orig_sum, final_sum)

        if PLOTS:
            for fov in opt.fov_manager.fovs:
                cnrs = fov.corners[1]
                plt.plot(cnrs[0], cnrs[1])

            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm(),
                       extent=(-opt.image_plane.hdu.header["NAXIS1"] / 2,
                               opt.image_plane.hdu.header["NAXIS1"] / 2,
                               -opt.image_plane.hdu.header["NAXIS2"] / 2,
                               opt.image_plane.hdu.header["NAXIS2"] / 2,))
            plt.colorbar()
            plt.show()

    def test_observe_works_for_many_sources_distributed(self, cmds, im_src):
        orig_sum = np.sum(im_src.fields[0].data)
        im_src.fields[0].data += 1
        im_src1 = deepcopy(im_src)
        im_src2 = deepcopy(im_src)
        im_src2.shift(7, 7)
        im_src3 = deepcopy(im_src)
        im_src3.shift(-10, 14)
        im_src4 = deepcopy(im_src)
        im_src4.shift(-4, -6)
        im_src5 = deepcopy(im_src)
        im_src5.shift(15, -15)
        multi_img = im_src1 + im_src2 + im_src3 + im_src4 + im_src5

        cmds["SIM_PIXEL_SCALE"] = 0.02
        opt = OpticalTrain(cmds)
        opt.observe(multi_img)

        final_sum = np.sum(opt.image_plane.image)
        print(orig_sum, final_sum)

        if PLOTS:
            for fov in opt.fov_manager.fovs:
                cnrs = fov.corners[1]
                plt.plot(cnrs[0], cnrs[1])

            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm(),
                       extent=(-opt.image_plane.hdu.header["NAXIS1"] / 2,
                               opt.image_plane.hdu.header["NAXIS1"] / 2,
                               -opt.image_plane.hdu.header["NAXIS2"] / 2,
                               opt.image_plane.hdu.header["NAXIS2"] / 2,))
            plt.colorbar()
            plt.show()


@pytest.mark.usefixtures("unity_cmds", "unity_src")
class TestReadout:
    def test_readout_works_when_source_observed(self, unity_cmds, unity_src):

        opt = OpticalTrain(unity_cmds)
        opt.observe(unity_src)
        hdu = opt.readout()

        if PLOTS:
            plt.subplot(221)
            plt.imshow(unity_src.fields[0].data)
            plt.colorbar()

            plt.subplot(222)
            plt.imshow(opt.image_plane.image)
            plt.colorbar()

            plt.subplot(223)
            plt.imshow(hdu[1].data)
            plt.colorbar()
            plt.show()

        src_average = np.average(unity_src.fields[0].data)
        assert np.average(hdu[1].data) == approx(src_average, rel=1e-3)
