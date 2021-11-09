import os
import pytest

import numpy as np

from scopesim.optics.fov2 import FieldOfView
from scopesim.optics.fov_manager2 import FOVManager

from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim.tests.mocks.py_objects import yaml_objects as yo
from scopesim.tests.mocks.py_objects import integr_spectroscopy_objects as iso

import matplotlib.pyplot as plt

PLOTS = False


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(preload_fovs=False), FOVManager)

    def test_initialises_with_list_of_effects(self):
        effects = eo._mvs_effects_list()
        assert isinstance(FOVManager(effects, preload_fovs=False), FOVManager)


class TestGenerateFovList:
    def test_returns_default_single_entry_fov_list_for_no_effects(self):
        fov_man = FOVManager(pixel_scale=1, plate_scale=1)
        fovs = fov_man.generate_fovs_list()
        assert len(fovs) == 1

    def test_returns_single_fov_for_mvs_system(self):
        effects = eo._mvs_effects_list()
        fov_man = FOVManager(effects=effects, pixel_scale=1, plate_scale=1)
        fovs = fov_man.generate_fovs_list()
        fov_volume = fovs[0].volume()

        for effect in fov_man.effects:
            print(effect)

        assert len(fovs) == 1
        assert fov_volume["xs"][0] == -1024 / 3600      # [deg] 2k detector / pixel_scale
        assert fov_volume["waves"][0] == 0.6            # [um] filter blue edge