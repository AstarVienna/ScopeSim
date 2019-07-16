import os
import pytest

from scopesim import rc
from scopesim.commands.user_commands2 import UserCommands

rc.__config__["!SIM.file.local_packages_path"] = "C:/Work/irdb/"


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
        cmd = UserCommands(use_instrument="MICADO")
        for yaml_dic in cmd.yaml_dicts:
            print(yaml_dic)
        print(cmd)


