import os
import shutil
import pytest

from scopesim import rc
from scopesim.commands.user_commands import UserCommands
from scopesim.server.database import download_package

LOCAL_PKGS_PATH = "./scopesim_pkg_dir_tmp/"
rc.__config__["!SIM.file.local_packages_path"] = os.path.abspath(LOCAL_PKGS_PATH)
# rc.__config__["!SIM.file.local_packages_path"] = "D:/Work/irdb/"


def setup_module():
    if not os.path.exists(LOCAL_PKGS_PATH):
        os.mkdir(LOCAL_PKGS_PATH)
    download_package("instruments/test_package.zip")


def teardown_module():
    shutil.rmtree(LOCAL_PKGS_PATH)


class TestInit:
    def test_initialise_with_nothing(self):
        assert isinstance(UserCommands(), UserCommands)

    def test_initialise_with_correct_keywords(self):
        cmd = UserCommands(yamls=["test_instrument.yaml"],
                           packages=["test_package"],
                           properties={"life": 42})
        assert isinstance(cmd, UserCommands)
        assert cmd["!OBS.life"] == 42
        assert cmd["!INST.pixel_scale"] == 0.5
        print(cmd)

    def test_initialise_with_keyword_use_instrument(self):
        cmd = UserCommands(use_instrument="test_package")
        assert isinstance(cmd, UserCommands)
        assert cmd["!INST.pixel_scale"] == 0.5
        print(cmd)


class TestUpdate:
    def test_updates_with_yaml_dict(self):
        yaml_input = {"alias": "TEL",
                      "properties": {"temperature": 8999}}
        cmd = UserCommands(use_instrument="test_package")
        cmd.update(yaml_input)
        assert cmd["!TEL.temperature"] < 9000

    def test_update_works_via_setitem(self):
        cmd = UserCommands(use_instrument="test_package")
        cmd["!TEL.gigawatts"] = 1.21
        assert cmd["!TEL.gigawatts"] == 1.21


class TestYamlDicts:
    def test_everything_original_is_in_the_yaml_dicts_list(self):
        cmd = UserCommands(use_instrument="test_package")
        assert cmd["!TEL.temperature"] > 9000
        assert len(cmd.yaml_dicts) > 0
        # for yaml_dic in cmd.yaml_dicts:
        #     print(yaml_dic)


class TestListLocalPackages:
    def test_all_packages_listed(self):
        from scopesim.commands import user_commands as uc2
        real_pkgs, ext_pkgs = uc2.list_local_packages(action="return")
        assert len(real_pkgs) > 0
