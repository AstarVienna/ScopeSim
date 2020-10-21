import os
import pytest
from pytest import approx

import numpy as np

from scopesim import rc
from scopesim.optics.optical_train import OpticalTrain
from scopesim.utils import find_file
from scopesim.commands import UserCommands

from scopesim.tests.mocks.py_objects.source_objects import _image_source, \
    _single_table_source

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


@pytest.fixture(scope="function")
def cmds():
    return UserCommands(yamls=["CMD_unity_cmds.yaml"])


@pytest.fixture(scope="function")
def non_unity_cmds():
    return UserCommands(yamls=["CMD_non_unity_cmds.yaml"])


@pytest.fixture(scope="function")
def tbl_src():
    return _single_table_source()


@pytest.fixture(scope="function")
def im_src():
    return _image_source()


class TestObserve:
    # The CMD_unity_cmds.yaml sets the background emission to 0
    @pytest.mark.usefixtures("cmds", "im_src", "tbl_src")
    def test_flux_is_conserved_for_no_bg_emission(self, cmds, tbl_src):
        opt = OpticalTrain(cmds)
        opt.observe(tbl_src)
        im = opt.image_planes[0].image
        bg_flux = np.pi / 4 * np.prod(im.shape)
        src_flux = tbl_src.photons_in_range(1, 2, 1)[0].value

        if PLOTS:
            implane = opt.image_planes[0]
            plt.imshow(implane.image.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

        # given a 1 um bandpass
        print(src_flux, bg_flux)
        area = opt.optics_manager.area.value  # u.m**2
        assert src_flux == approx(1)          # u.Unit("ph s-1")
        assert np.sum(im) == approx(src_flux * area, rel=2e-3)

    @pytest.mark.usefixtures("non_unity_cmds", "tbl_src")
    def test_flux_is_conserved_for_yes_bg_emission(self, non_unity_cmds, tbl_src):
        """
        # originally all has 1 count. atmo TC=0.9, mirror TC=0.5. Hence:
        # 0.5 counts from atmo, 1 count from mirror, 0.45 count from source
        """
        opt = OpticalTrain(non_unity_cmds)
        opt.observe(tbl_src)
        im = opt.image_planes[0].image
        src_flux = tbl_src.photons_in_range(1, 2, 1)[0].value

        if PLOTS:
            implane = opt.image_planes[0]
            plt.imshow(implane.image.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

        # given a 1 um bandpass
        assert src_flux == approx(1)          # u.Unit("ph s-1")
        assert np.max(im) - np.median(im) == approx(0.45, rel=1e-2)
        assert np.median(im) == approx(1.5, abs=1e-2)
