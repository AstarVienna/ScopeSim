# integration test using everything and the HAWKI package
import pytest
from pytest import approx
import os
import shutil

import numpy as np
from astropy import units as u

import scopesim
import scopesim.source.source_templates
from scopesim import rc

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

# pytest_plugins = ['pytest_profiling']
if rc.__config__["!SIM.tests.ignore_integration_tests"]:
    pytestmark = pytest.mark.skip("Ignoring HAWKI integration tests")

rc.__config__["!SIM.file.local_packages_path"] = "./scopesim_pkg_dir_tmp/hawki_tmp/"

PKGS = {"Paranal": "locations/Paranal.zip",
        "VLT": "telescopes/VLT.zip",
        "HAWKI": "instruments/HAWKI.zip"}

CLEAN_UP = True
PLOTS = False


def setup_module():
    rc.__config__["!SIM.file.use_cached_downloads"] = False
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
        assert rc.__currsys__["!TEL.area"].value == approx(52.81, rel=1e-3)
        assert rc.__currsys__["!TEL.etendue"].value == approx(0.5934, rel=1e-3)
        assert rc.__currsys__["!INST.pixel_scale"] == approx(0.106, rel=1e-3)

        # test that OpticalTrain builds properly
        assert isinstance(opt, scopesim.OpticalTrain)

        # test that we have a system throughput
        wave = np.linspace(0.7, 2.5, 181) * u.um
        tc = opt.optics_manager.surfaces_table.throughput
        ec = opt.optics_manager.surfaces_table.emission
        # ..todo:: something super wierd is going on here when running pytest in the top directory
        assert 0.55 < np.max(tc(wave)) < 0.7

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


class TestObserveOpticalTrain:
    def test_background_is_similar_to_online_etc(self):
        cmd = scopesim.UserCommands(use_instrument="HAWKI")
        opt = scopesim.OpticalTrain(cmd)
        src = scopesim.source.source_templates.empty_sky()

        # ETC gives 2700 e-/DIT for a 1s DET at airmass=1.2, pwv=2.5
        opt.observe(src)
        assert np.average(opt.image_planes[0].data) == approx(2700, rel=0.2)

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
