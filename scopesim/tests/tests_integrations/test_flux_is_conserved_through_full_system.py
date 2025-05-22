
import pytest
from pytest import approx
from unittest.mock import patch

import numpy as np

from scopesim.optics.optical_train import OpticalTrain
from scopesim.commands import UserCommands

from scopesim.tests.mocks.py_objects.source_objects import (
    _image_source, _single_table_source, _table_source_overlapping)

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


PLOTS = False


@pytest.fixture(scope="function")
def cmds(mock_path, mock_path_yamls):
    with patch("scopesim.rc.__search_path__", [mock_path, mock_path_yamls]):
        return UserCommands(yamls=["CMD_unity_cmds.yaml"])


@pytest.fixture(scope="function")
def non_unity_cmds(mock_path, mock_path_yamls):
    with patch("scopesim.rc.__search_path__", [mock_path, mock_path_yamls]):
        return UserCommands(yamls=["CMD_non_unity_cmds.yaml"])


@pytest.fixture(scope="function")
def tbl_src():
    return _single_table_source()


@pytest.fixture(scope="function")
def im_src():
    return _image_source()


@pytest.mark.usefixtures("protect_currsys", "patch_mock_path")
class TestObserve:
    # The CMD_unity_cmds.yaml sets the background emission to 0
    def test_flux_is_conserved_for_no_bg_emission(self, cmds, tbl_src):
        opt = OpticalTrain(cmds)
        opt.observe(tbl_src)
        im = opt.image_planes[0].image
        bg_flux = np.pi / 4 * np.prod(im.shape)
        src_flux = tbl_src.photons_in_range(1, 2, 1)[0].value

        if PLOTS:
            implane = opt.image_planes[0]
            plt.imshow(implane.image.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

        # given a 1 um bandpass
        print(src_flux, bg_flux)
        area = opt.optics_manager.area.value  # u.m**2
        assert src_flux == approx(1)          # u.Unit("ph s-1")
        assert np.sum(im) == approx(src_flux * area, rel=2e-3)

    def test_flux_is_conserved_for_yes_bg_emission(self, non_unity_cmds, tbl_src):
        """
        # originally all has 1 count. atmo TC=0.9, mirror TC=0.5. Hence:
        # 0.5 counts from atmo, 1 count from mirror, 0.45 count from source
        """
        opt = OpticalTrain(non_unity_cmds)
        opt.observe(tbl_src)
        im = opt.image_planes[0].image
        src_flux = tbl_src.photons_in_range(1, 2, 1)[0].value

        if PLOTS:
            implane = opt.image_planes[0]
            plt.imshow(implane.image.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

        # given a 1 um bandpass
        assert src_flux == approx(1)          # u.Unit("ph s-1")
        assert np.sum(im - np.median(im)) == approx(0.45, rel=1e-2)
        assert np.median(im) == approx(1.5, abs=1e-2)


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestStackedStars:
    """Test whether stars can be stacked."""

    def test_stacked_stars(self):
        """Test whether stars can be stacked.

        Four stacked faint stars should have the same magnitude as one bright star.
        """
        stars = _table_source_overlapping()
        temp_celsius = 5.0
        dit = 10  # s
        ndit = 1
        s_filter = "H"
        cmd = UserCommands(use_instrument="basic_instrument")

        cmd.update(properties={
            "!OBS.ndit": ndit,
            "!OBS.dit": dit,
            "!OBS.airmass": 1.0,
            "!ATMO.background.filter_name": s_filter,
            "!ATMO.temperature": temp_celsius,
            "!TEL.temperature": temp_celsius,
            "!DET.width": 512,  # pixel
            "!DET.height": 512,
        })
        micado = OpticalTrain(cmd)

        # disabling effects
        micado["dark_current"].include = False
        micado["shot_noise"].include = False
        micado["detector_linearity"].include = False
        micado["exposure_action"].include = False
        micado['readout_noise'].include = False
        micado["source_fits_keywords"].include = False
        micado["effects_fits_keywords"].include = False
        micado["config_fits_keywords"].include = False
        micado["extra_fits_keywords"].include = False
        micado["telescope_reflection"].include = False
        micado["qe_curve"].include = False
        micado["static_surfaces"].include = False
        micado.observe(stars)
        hdus_h = micado.readout()

        im_h = hdus_h[0][1].data
        dx, dy = 255, 255

        if PLOTS:
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 12))
            axes.imshow(np.sqrt(im_h), norm=LogNorm(), cmap="inferno")
            axes.set_title('H-band MORFEO "PSF Generic" ')
            plt.show()

        quadrants = im_h[:dx, :dy], im_h[dx:, :dy], im_h[:dx, dy:], im_h[dx:, dy:]
        the_sums = [q.sum() for q in quadrants]
        flux_one_star, empty1, empty2, flux_four_stacked_stars = the_sums

        # Check whether the stars are equal.
        assert flux_one_star == pytest.approx(flux_four_stacked_stars, rel=0.05)

        # Check whether the empty skies are equal.
        assert empty1 == pytest.approx(empty2, rel=0.05)

        # Check whether the star is brighter than the sky.
        assert flux_one_star > empty1 * 2
