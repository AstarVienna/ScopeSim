# integration test using everything and the HAWKI package
import pytest
from pytest import approx
import os
import shutil

import numpy as np
from astropy import units as u
from astropy.io import ascii

import scopesim
import scopesim.source.source_templates
from scopesim.tests.mocks.py_objects.source_objects import _single_table_source
from scopesim import rc

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

# pytest_plugins = ['pytest_profiling']
if rc.__config__["!SIM.tests.run_integration_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring HAWKI integration tests")


PKGS = {"Paranal": "locations/Paranal.zip",
        "VLT": "telescopes/VLT.zip",
        "HAWKI": "instruments/HAWKI.zip"}

CLEAN_UP = True
PLOTS = False


def setup_module():
    rc.__config__["!SIM.file.use_cached_downloads"] = False
    rc_local_path = "./TEMP_HAWKI/"
    rc.__config__["!SIM.file.local_packages_path"] = rc_local_path

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
        print("irdb" not in rc_local_path)
        for pkg_name in PKGS:
            assert os.path.isdir(os.path.join(rc_local_path, pkg_name))


class TestLoadUserCommands:
    def test_user_commands_loads_without_throwing_errors(self, capsys):
        cmd = scopesim.UserCommands(use_instrument="HAWKI")
        assert isinstance(cmd, scopesim.UserCommands)
        for key in ["SIM", "OBS", "ATMO", "TEL", "INST", "DET"]:
            assert key in cmd and len(cmd[key]) > 0

        stdout = capsys.readouterr()
        assert len(stdout.out) == 0


class TestMakeOpticalTrain:
    def test_works_seamlessly_for_hawki_package(self, capsys):
        cmd = scopesim.UserCommands(use_instrument="HAWKI")
        opt = scopesim.OpticalTrain(cmd)

        # test that the major values have been updated in rc.__currsys__
        assert rc.__currsys__["!TEL.area"].value == approx(52.02, rel=1e-3)
        assert rc.__currsys__["!TEL.etendue"].value == approx(0.58455, rel=1e-3)
        assert rc.__currsys__["!INST.pixel_scale"] == approx(0.106, rel=1e-3)

        # test that OpticalTrain builds properly
        assert isinstance(opt, scopesim.OpticalTrain)

        # test that we have a system throughput
        wave = np.linspace(0.7, 2.5, 181) * u.um
        tc = opt.optics_manager.surfaces_table.throughput
        ec = opt.optics_manager.surfaces_table.emission
        # ..todo:: something super wierd is going on here when running pytest in the top directory
        assert 0.5 < np.max(tc(wave)).value < 0.55

        if PLOTS:
            plt.plot(wave, tc(wave))
            plt.show()

        if PLOTS:
            plt.plot(wave, ec(wave))
            plt.show()

        # test that we have the correct number of FOVs for Ks band
        assert len(opt.fov_manager.fovs) == 9

        if PLOTS:
            fovs = opt.fov_manager.fovs
            from scopesim.optics.image_plane_utils import calc_footprint
            plt.subplot(121)
            for fov in fovs:
                x, y = calc_footprint(fov.hdu.header)
                plt.fill(x*3600, y*3600, alpha=0.1, c="b")
                plt.title("Sky plane")
                plt.xlabel("[arcsec]")

            plt.subplot(122)
            for fov in fovs:
                x, y = calc_footprint(fov.hdu.header, "D")
                plt.fill(x, y)
                plt.title("Detector focal plane")
                plt.xlabel("[mm]")

            plt.show()

        # test that the ImagePlane is large enough
        assert opt.image_planes[0].header["NAXIS1"] > 4200
        assert opt.image_planes[0].header["NAXIS2"] > 4200
        assert np.all(opt.image_planes[0].data == 0)

        # test assert there are 4 detectors, each 2048x2048 pixels
        hdu = opt.readout()[0]
        assert len(opt.detector_arrays[0].detectors) == 4
        for detector in opt.detector_arrays[0].detectors:
            assert detector.hdu.header["NAXIS1"] == 2048
            assert detector.hdu.header["NAXIS2"] == 2048

        if PLOTS:
            for i in range(1, 5):
                plt.subplot(2, 2, i)
                plt.imshow(hdu[i].data)
            plt.show()

        dit = rc.__currsys__["!OBS.dit"]
        ndit = rc.__currsys__["!OBS.ndit"]
        assert np.average(hdu[1].data) == approx(ndit * dit * 0.1, abs=0.5)

    def test_system_transmission_is_similar_to_eso_etc(self):
        """
        A ~20% discrepency between the ESO and ScopeSim system throughputs
        """

        for filt_name in ["Y", "J", "H", "Ks", "BrGamma", "CH4"]:

            cmd = scopesim.UserCommands(use_instrument="HAWKI")
            cmd["!OBS.filter_name"] = filt_name
            opt = scopesim.OpticalTrain(cmd)
            opt["paranal_atmo_default_ter_curve"].include = False

            src = _single_table_source(n=1000)
            opt.observe(src)

            if PLOTS:
                fname = "hawki_eso_etc/TER_system_{}.dat".format(filt_name)
                dname = os.path.dirname(__file__)
                etc_tbl = ascii.read(os.path.join(dname, fname))
                etc_wave = etc_tbl["wavelength"] * 1e-3 * u.um
                etc_thru = etc_tbl["transmission"] * 1e-2
                # plt.plot(etc_wave, etc_thru, c="b", label="ESO/ETC")

                flux_init = src.spectra[0](etc_wave)
                flux_final = opt._last_source.spectra[0](etc_wave)
                ss_thru = flux_final / flux_init
                # plt.plot(etc_wave, ss_thru, c="r", label="ScopeSim/HAWKI")

                plt.plot(etc_wave, ss_thru / etc_thru - 1)

                plt.ylim(0, 0.5)
                plt.show()


class TestObserveOpticalTrain:
    def test_background_is_similar_to_online_etc(self):
        """
        Based on mocks/photometry/check_photometry.py

        Sky BG from my HAWKI archive data (ph s-1 pixel-1)
        J 170
        H 958
        Ks 1204
        BrG 213

        K filter
        --------
        Skycalc BG ph flux for K: 0.08654 ph / (cm2 s arcsec2)
        -> HAWKI sky BG = 510 ph / s / pixel
                        = 0.08654252 * (410**2 * np.pi) * (0.106**2)

        ETC gives 2550 e-/DIT/pixel for a 1s DET at airmass=1.0, pwv=2.5
        (HAWKI archive average 1204 e-/s-1/pixel)
        Remaining must come from VLT or entrance window

        ScopeSim Flux contributors to final BG count
        - all : 1360
        - minus "paranal_atmo_default_ter_curve" : 850
        - minus "vlt_mirror_list" : 780
        - minus entrance window from "hawki_mirror_list" : 0

        H filter
        --------
        Skycalc BG ph flux for K: 0.39862 ph / (cm2 s arcsec2)
        -> HAWKI sky BG = 2365 ph / s / pixel

        ETC gives 1625 e-/DIT/pixel for a 1s DET at airmass=1.0, pwv=2.5
        (HAWKI archive average 958 e-/s-1/pixel)
        Remaining must come from VLT or entrance window

        ScopeSim Flux contributors to final BG count
        - all : 2370
        - minus "paranal_atmo_default_ter_curve" : 2
        - minus "vlt_mirror_list" : 1.8
        - minus entrance window from "hawki_mirror_list" : 0

        J filter
        --------
        Skycalc BG ph flux for K: 0.39862 ph / (cm2 s arcsec2)
        -> HAWKI sky BG = 444 ph / s / pixel

        ETC gives 225 e-/DIT/pixel for a 1s DET at airmass=1.0, pwv=2.5
        (HAWKI archive average 170 e-/s-1/pixel)
        Remaining must come from VLT or entrance window

        ScopeSim Flux contributors to final BG count
        - all : 290
        - minus "paranal_atmo_default_ter_curve" : 0
        - minus "vlt_mirror_list" : 0
        - minus entrance window from "hawki_mirror_list" : 0

        Given the variability of the backgrounds, ScopeSim is doing a pretty
        good job of making the background contributions


        """
        cmd = scopesim.UserCommands(use_instrument="HAWKI")
        cmd["!OBS.filter_name"] = "H"
        opt = scopesim.OpticalTrain(cmd)
        opt["paranal_atmo_default_ter_curve"].include = True
        opt["vlt_mirror_list"].include = True
        opt["hawki_mirror_list"].include = True

        src = scopesim.source.source_templates.empty_sky()
        opt.observe(src)

        if PLOTS:
            wave = np.arange(0.7, 2.5, 0.001) * u.um
            specs = opt._last_source.spectra
            for i in range(len(specs)):
                flux = specs[i](wave)
                plt.plot(wave, flux)
            plt.show()

        assert np.average(opt.image_planes[0].data) == approx(2400, rel=0.2)

    def test_actually_produces_stars(self):
        cmd = scopesim.UserCommands(use_instrument="HAWKI",
                                    properties={"!OBS.dit": 360,
                                                "!OBS.ndit": 10})
        cmd.ignore_effects += ["detector_linearity"]

        opt = scopesim.OpticalTrain(cmd)
        src = scopesim.source.source_templates.star_field(10000, 5, 15, 440)

        # ETC gives 2700 e-/DIT for a 1s DET at airmass=1.2, pwv=2.5
        # background should therefore be ~ 8.300.000
        opt.observe(src)
        hdu = opt.readout()[0]

        implane_av = np.average(opt.image_planes[0].data)
        hdu_av = np.average([hdu[i].data for i in range(1, 5)])
        exptime = cmd["!OBS.ndit"] * cmd["!OBS.dit"]

        assert hdu_av == approx(implane_av * exptime, rel=0.01)

        if PLOTS:
            plt.subplot(1, 2, 1)
            plt.imshow(opt.image_planes[0].image[128:2048, 128:2048].T,
                       norm=LogNorm())
            plt.colorbar()

            plt.subplot(1, 2, 2)
            plt.imshow(hdu[1].data[128:2048, 128:2048].T, norm=LogNorm(),
                       vmax=3e7)
            plt.colorbar()

            plt.show()
