import os
import pytest
from pytest import approx

import numpy as np

import scopesim as sim
from scopesim.optics.optical_train import OpticalTrain
from scopesim.utils import find_file
from scopesim.commands.user_commands2 import UserCommands

from scopesim.tests.mocks.py_objects.source_objects import _image_source, \
    _single_table_source

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))

sim.rc.__search_path__ += [FILES_PATH, YAMLS_PATH]


def _basic_cmds():
    return UserCommands(filename=find_file("CMD_unity_cmds.config"))


@pytest.fixture(scope="function")
def cmds():
    return _basic_cmds()


@pytest.fixture(scope="function")
def tbl_src():
    return _single_table_source()


@pytest.fixture(scope="function")
def im_src():
    return _image_source()


@pytest.mark.usefixtures("cmds", "im_src", "tbl_src")
class TestObserve:
    def test_flux_is_conserved_and_emission_level_correct(self, cmds, tbl_src):
        opt = OpticalTrain(cmds)
        opt.observe(tbl_src)
        im = opt.image_plane.image
        bg_flux = np.pi / 4 * np.prod(im.shape)
        src_flux = tbl_src.photons_in_range(1, 2, 1)[0].value

        if PLOTS is False:
            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

        assert src_flux == approx(1)          # u.Unit("ph s-1")
        assert np.sum(im) == approx(src_flux + bg_flux, rel=2e-3)
        print(src_flux, bg_flux)

        # given a 1 um bandpass
