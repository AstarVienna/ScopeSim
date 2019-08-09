# integration test using everything and the MICADO package
import pytest
import os
import shutil

import scopesim
from scopesim import rc

CLEAN_UP = True
rc.__config__["!SIM.file.local_packages_path"] = "./scopesim_pkg_dir_tmp/"
# rc.__config__["!SIM.file.local_packages_path"] = "C:/Work/irdb/"

PKGS = {"Armazones": "locations/Armazones.zip",
        "ELT": "telescopes/ELT.zip",
        "MAORY": "instruments/MAORY.zip",
        "MICADO": "instruments/MICADO.zip"}


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
        for pkg_name in PKGS:
            rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
            assert os.path.isdir(os.path.join(rc_local_path, pkg_name))
        print("irdb" not in rc_local_path)


class TestLoadUserCommands:
    def test_user_commands_loads_with_throwing_errors(self, capsys):
        cmd = scopesim.UserCommands(use_instrument="MICADO")
        assert isinstance(cmd, scopesim.UserCommands)
        for key in ["SIM", "OBS", "ATMO", "TEL", "INST", "RO", "INST", "DET"]:
            assert key in cmd and len(cmd[key]) > 0

        stdout = capsys.readouterr()
        assert len(stdout.out) == 0


class TestMakeOpticalTrain:
    def test_works_seamlessly_for_micado_package(self, capsys):
        cmd = scopesim.UserCommands(use_instrument="MICADO")
        opt = scopesim.OpticalTrain(cmd)
        assert isinstance(opt, scopesim.OpticalTrain)

        stdout = capsys.readouterr()
        assert len(stdout.out) == 0
