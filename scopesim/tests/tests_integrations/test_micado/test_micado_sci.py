# integration test using everything and the MICADO package
import pytest
from pytest import approx
import os
import shutil

import numpy as np
from astropy import units as u
from astropy.io import fits

import scopesim as sim
from scopesim import rc

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

CLEAN_UP = True
PLOTS = False

rc.__config__["!SIM.file.local_packages_path"] = "C:/Work/irdb/"


class TestObserve:
    def test_star_field_with_mcao_4mas(self):
        cmd = sim.UserCommands(use_instrument="MICADO_Sci",
                               set_modes=["SCAO", "4mas"])
        opt = sim.OpticalTrain(cmd)

        src = sim.source.source_templates.star_field(100, 20, 30, 3, use_grid=True)
        opt.observe(src)

        plt.imshow(opt.image_planes[0].image, norm=LogNorm())
        plt.show()

    def test_star_field_with_spec(self):
        cmd = sim.UserCommands(use_instrument="MICADO_Sci",
                               set_modes=["MCAO", "SPEC"])
        opt = sim.OpticalTrain(cmd)

        src = sim.source.source_templates.star_field(100, 20, 30, 3, use_grid=True)
        opt.observe(src)

        plt.imshow(opt.image_planes[0].image, norm=LogNorm())
        plt.show()
