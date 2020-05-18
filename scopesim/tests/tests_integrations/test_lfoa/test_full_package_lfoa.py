# integration test using everything and the LFAO package
import pytest
from pytest import approx
import os
import shutil

import numpy as np
from astropy import units as u
from astropy.io import fits

import scopesim
from scopesim import rc

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


# if rc.__config__["!SIM.tests.run_integration_tests"] is False:
#     pytestmark = pytest.mark.skip("Ignoring MICADO integration tests")

rc.__config__["!SIM.file.local_packages_path"] = "./lfoa_temp/"
rc.__config__["!SIM.file.use_cached_downloads"] = False

PKGS = {"LFOA": "telescopes/LFOA.zip"}

CLEAN_UP = True
PLOTS = False


def setup_module():
    rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
    if not os.path.exists(rc_local_path):
        os.mkdir(rc_local_path)
        rc.__config__["!SIM.file.local_packages_path"] = os.path.abspath(
            rc_local_path)

    for pkg_name in PKGS:
        if not os.path.isdir(os.path.join(rc_local_path, pkg_name)) and \
                "irdb" not in rc_local_path:
            scopesim.download_package(PKGS[pkg_name])


def teardown_module():
    rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
    if CLEAN_UP and "irdb" not in rc_local_path:
        shutil.rmtree(rc.__config__["!SIM.file.local_packages_path"])


class TestInit:
    def test_all_packages_are_available(self):
        rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
        for pkg_name in PKGS:
            assert os.path.isdir(os.path.join(rc_local_path, pkg_name))
        print("irdb" not in rc_local_path)


class TestLoadUserCommands:
    def test_user_commands_loads_without_throwing_errors(self, capsys):
        cmd = scopesim.UserCommands(use_instrument="LFOA")
        assert isinstance(cmd, scopesim.UserCommands)

        stdout = capsys.readouterr()
        assert len(stdout.out) == 0


class TestMakeOpticalTrain:
    def test_load_lfao(self):
        cmd = scopesim.UserCommands(use_instrument="LFOA",
                                    properties={"!OBS.filter_name": "OIII",
                                                "!OBS.dit": 0})
        opt = scopesim.OpticalTrain(cmd)
        assert isinstance(opt, scopesim.OpticalTrain)

        src = scopesim.source.source_templates.star_field(5000, 10, 20, 700)
        # src = scopesim.source.source_templates.empty_sky()
        opt.observe(src)
        hdu_list = opt.readout()[0]

        assert isinstance(hdu_list, fits.HDUList)

        plt.imshow(hdu_list[1].data, norm=LogNorm())
        plt.show()
