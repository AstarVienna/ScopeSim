import pytest
from pytest import raises
import os

import scopesim as sim
from scopesim.source import source_templates as st


inst_pkgs = os.path.join(os.path.dirname(__file__), "../mocks")
sim.rc.__currsys__["!SIM.file.local_packages_path"] = inst_pkgs

class TestLoadsUserCommands:
    def test_loads(self):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        assert isinstance(cmd, sim.UserCommands)
        assert cmd["!INST.pixel_scale"] == 0.2

class TestLoadsOpticalTrain:
    def test_loads(self):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        opt = sim.OpticalTrain(cmd)
        assert isinstance(opt, sim.OpticalTrain)
        assert opt["#slit_wheel.current_slit"] == "!OBS.slit_name"
        assert opt["#slit_wheel.current_slit!"] == "narrow"
