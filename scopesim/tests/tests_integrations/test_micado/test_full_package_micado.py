# integration test using everything and the MICADO package
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


if rc.__config__["!SIM.tests.ignore_integration_tests"]:
    pytestmark = pytest.mark.skip("Ignoring MICADO integration tests")

rc.__config__["!SIM.file.local_packages_path"] = "./scopesim_pkg_dir_tmp/"


PKGS = {"Armazones": "locations/Armazones.zip",
        "ELT": "telescopes/ELT.zip",
        "MAORY": "instruments/MAORY.zip",
        "MICADO": "instruments/MICADO.zip"}

CLEAN_UP = True
PLOTS = False


def setup_module():
    rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
    if not os.path.exists(rc_local_path) or CLEAN_UP:
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
        cmd = scopesim.UserCommands(use_instrument="MICADO")
        assert isinstance(cmd, scopesim.UserCommands)
        for key in ["SIM", "OBS", "ATMO", "TEL", "INST", "DET"]:
            assert key in cmd and len(cmd[key]) > 0

        stdout = capsys.readouterr()
        assert len(stdout.out) == 0

    def test_user_commands_loads_mode_files(self):
        cmd = scopesim.UserCommands(use_instrument="MICADO")
        yaml_names = [yd["name"] for yd in cmd.yaml_dicts]
        print(yaml_names)

        assert "MICADO_IMG_LR" in yaml_names

    def test_user_commands_can_change_modes(self):
        cmd = scopesim.UserCommands(use_instrument="MICADO")
        cmd.set_mode("mcao_spec")
        assert "MAORY" in [yd["name"] for yd in cmd.yaml_dicts]
        assert "MICADO_SPEC" in [yd["name"] for yd in cmd.yaml_dicts]

    def test_user_commands_can_change_modes_via_init(self):
        cmd = scopesim.UserCommands(use_instrument="MICADO",
                                    set_mode="mcao_spec")
        assert "MAORY" in [yd["name"] for yd in cmd.yaml_dicts]
        assert "MICADO_SPEC" in [yd["name"] for yd in cmd.yaml_dicts]


class TestMakeOpticalTrain:
    def test_works_seamlessly_for_micado_package(self, capsys):

        src = scopesim.source.source_templates.star_field(10000, 10, 25, 60)
        src = scopesim.source.source_templates.empty_sky()
        cmd = scopesim.UserCommands(use_instrument="MICADO")
        opt = scopesim.OpticalTrain(cmd)
        assert isinstance(opt, scopesim.OpticalTrain)

        # ..todo:: add source object here!
        # opt.observe(src)
        hdu = opt.readout()[0]
        assert isinstance(hdu, fits.HDUList)
