import pytest
from pytest import raises
import os
from time import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u

import scopesim as sim
from scopesim.source import source_templates as st

PLOTS = False

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


class TestObserveImagingMode:
    def test_runs(self):
        src = st.star(flux=9)
        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["imaging"])
        opt = sim.OpticalTrain(cmd)
        opt.observe(src)
        hdul = opt.readout()[0]

        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(opt.image_planes[0].data, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im, norm=LogNorm())
            plt.show()

        assert det_im.max() > 1e5
        assert det_im[505:520, 505:520].sum() > 3e6


class TestObserveSpectroscopyMode:
    """
    Test the number of spots along the three spectral traces.

    Spots are places at 0.05um intervals, offset by 0.025um

    3 traces vertically covering the regions:
    - J: 1.0, 1.35     dwave = 0.35 --> R~2800      -> 7 spots
    - H: 1.3, 1.8      dwave = 0.5  --> R~2000      -> 10 spots
    - K: 1.75, 2.5     dwave = 0.75 --> R~1300      -> 15 spots
    """
    def test_runs(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.zeros(len(wave))
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=[0], y=[0], ref=[0], weight=[1e-3])

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["spectroscopy"])
        opt = sim.OpticalTrain(cmd)
        for effect_name in ["shot_noise", "dark_current", "readout_noise",
                            "atmospheric_radiometry", "source_fits_keywords",
                            "effects_fits_keywords", "config_fits_keywords"]:
            opt[effect_name].include = False

        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im)
            plt.show()

        xs = [(175, 200), (500, 525), (825, 850)]
        dlams = np.array([0.35, 0.5, 0.75])
        n_spots = dlams / 0.05

        spot_flux = 28000       # due to psf flux losses in narrow slit (0.5")
        for i in range(3):
            x0, x1 = xs[i]
            trace_flux = det_im[:, x0:x1].sum()     # sum along a trace
            assert round(trace_flux / spot_flux) == round(n_spots[i])


class TestObserveIfuMode:
    def test_runs(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.zeros(len(wave))
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        x = [-4, -2, 0, 2, 4] * 5
        y = [y0 for y0 in range(-4, 5, 2) for i in range(5)]
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=x, y=y, ref=[0]*len(x), weight=[1e-3]*len(x))

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["ifu"])
        opt = sim.OpticalTrain(cmd)
        for effect_name in ["shot_noise", "dark_current", "readout_noise",
                            "atmospheric_radiometry", "source_fits_keywords",
                            "effects_fits_keywords", "config_fits_keywords"]:
            opt[effect_name].include = False

        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im)
            plt.show()

        xs = [(157, 226), (317, 386), (476, 545), (640, 706), (797, 866)]
        spot_flux = 100000       # due to psf flux in large IFU slices (2")
        for i in range(5):
            x0, x1 = xs[i]
            trace_flux = det_im[:, x0:x1].sum()     # sum along a trace
            assert round(trace_flux / spot_flux) == 15 * 5


    def test_random_star_field(self):
        src = sim.source.source_templates.star_field(n=100, mmin=8, mmax=18, width=10)

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["ifu"])
        opt = sim.OpticalTrain(cmd)
        for effect_name in ["source_fits_keywords", "effects_fits_keywords",
                            "config_fits_keywords"]:
            opt[effect_name].include = False

        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im, norm=LogNorm())
            plt.show()


class TestFitsHeader:
    def test_source_keywords_in_header(self):
        src = st.star()
        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["imaging"])
        opt = sim.OpticalTrain(cmd)
        opt.observe(src)
        hdul = opt.readout()[0]
        hdr = hdul[0].header

        assert hdr["SIM SRC0 object"] == 'star'
        assert hdr["SIM EFF14 class"] == 'SourceDescriptionFitsKeywords'
        assert hdr["SIM CONFIG OBS filter_name"] == 'J'
        assert hdr["ESO ATM SEEING"] == sim.utils.from_currsys("!OBS.psf_fwhm")
