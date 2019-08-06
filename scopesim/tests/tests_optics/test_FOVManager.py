import os
import pytest
from pytest import approx

import numpy as np

import scopesim as sim
from scopesim import rc
from scopesim.optics.fov_manager import FOVManager
from scopesim.optics import fov_manager as fov_mgr
from scopesim.optics.image_plane import ImagePlane

from scopesim.tests.mocks.py_objects.effects_objects import _mvs_effects_list, \
    _atmospheric_dispersion, _filter_tophat_curve, _const_psf, _ncpa_psf
from scopesim.tests.mocks.py_objects.yaml_objects import \
     _usr_cmds_min_viable_scope

import matplotlib.pyplot as plt

PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__ += [NEW_PATH]


@pytest.fixture(scope="function")
def mvs_effects_list():
    return _mvs_effects_list()


@pytest.fixture(scope="function")
def mvs_usr_cmds():
    return _usr_cmds_min_viable_scope()


@pytest.mark.usefixtures("mvs_effects_list")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(), FOVManager)

    def test_initialises_with_list_of_effects(self, mvs_effects_list):
        assert isinstance(FOVManager(mvs_effects_list), FOVManager)


@pytest.mark.usefixtures("mvs_effects_list", "mvs_usr_cmds")
class TestGenerateFovs:
    def test_returns_the_desired_number_of_fovs(self, mvs_effects_list,
                                                mvs_usr_cmds):
        for yaml_dic in mvs_usr_cmds:
            rc.__currsys__.cmds.update(yaml_dic)
        fov_man = FOVManager(mvs_effects_list)
        fovs = fov_man.generate_fovs_list()

        implane_size = mvs_effects_list[-2].image_plane_header["NAXIS1"]
        chunk_size = sim.rc.__currsys__["!SIM.computing.chunk_size"]
        wave_layers = 1

        assert len(fovs) == np.round(implane_size / chunk_size)**2 * wave_layers

        if PLOTS:
            plt.subplot(121)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.hdu.header)
                plt.plot(x, y)

            plt.subplot(122)
            for fov in fovs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(fov.hdu.header, "D")
                plt.plot(x, y)

            plt.show()

    def test_fovs_dont_overlap_on_canvas(self, mvs_effects_list, mvs_usr_cmds):

        for yaml_dic in mvs_usr_cmds:
            rc.__currsys__.cmds.update(yaml_dic)

        implane = ImagePlane(mvs_effects_list[-2].image_plane_header)
        fov_man = FOVManager(mvs_effects_list)
        fovs = fov_man.generate_fovs_list()
        for fov in fovs:
            fov.hdu.data = np.ones((fov.header["NAXIS1"], fov.header["NAXIS2"]))
            implane.add(fov.hdu, wcs_suffix="D")

        if PLOTS:
            plt.imshow(implane.image.T, origin="lower")
            plt.colorbar()
            plt.show()

        # .. todo: Fix this - find out why!!!
        # assert np.all(implane.image == 1)


class TestGet3DShifts:
    def test_returns_zeros_when_no_shift3d_effects_passed(self):
        shifts = fov_mgr.get_3d_shifts([], wave_min=0.7, wave_max=3,
                                       pixel_scale=0.004, sub_pixel_fraction=1)
        assert np.all(shifts["x_shifts"] == 0)
        assert np.all(shifts["y_shifts"] == 0)

    def test_returns_almost_zero_for_zenith_atmospheric_dispersion(self):
        ad_zenith = _atmospheric_dispersion(airmass=1.)
        shifts = fov_mgr.get_3d_shifts([ad_zenith],
                                       pixel_scale=0.004, sub_pixel_fraction=1)
        assert np.all(shifts["x_shifts"] == 0)
        assert np.all(shifts["y_shifts"] == 0)

    def test_returns_non_zero_entries_for_off_zenith(self):
        ad_am_1_14 = _atmospheric_dispersion()
        shifts = fov_mgr.get_3d_shifts([ad_am_1_14],
                                       pixel_scale=0.004, sub_pixel_fraction=1)
        assert np.all(shifts["x_shifts"] == 0)
        assert np.interp(0, shifts["y_shifts"][::-1],
                         shifts["wavelengths"][::-1]) == pytest.approx(1.5)

    def test_shift_transferred_to_x_axis_for_90_deg_pupil_angle(self):
        ad_am_1_05 = _atmospheric_dispersion(airmass=1.05, pupil_angle=90)
        shifts = fov_mgr.get_3d_shifts([ad_am_1_05],
                                       pixel_scale=0.004, sub_pixel_fraction=1)
        assert np.interp(0, shifts["x_shifts"][::-1],
                         shifts["wavelengths"][::-1]) == pytest.approx(1.5)
        assert np.all(shifts["y_shifts"] == 0)

    def test_shifts_cancel_out_when_equal_and_opposite(self):
        ad_am_1_05 = _atmospheric_dispersion(airmass=1.05)
        ad_am_neg_1_04 = _atmospheric_dispersion(airmass=-1.05)
        shifts = fov_mgr.get_3d_shifts([ad_am_1_05, ad_am_neg_1_04],
                                       pixel_scale=0.004, sub_pixel_fraction=1)

        assert len(shifts["y_shifts"]) == 2
        assert np.min(shifts["wavelengths"]) == ad_am_1_05.meta["wave_min"]
        assert np.max(shifts["wavelengths"]) == ad_am_1_05.meta["wave_max"]

    @pytest.mark.parametrize("sub_pix_frac", [1, 0.3, 0.1])
    def test_combined_shifts_reduced_to_usable_number(self, sub_pix_frac):
        ad_am_1_05 = _atmospheric_dispersion(airmass=1.05)
        ad_am_neg_1_04 = _atmospheric_dispersion(airmass=-1.04)
        shifts = fov_mgr.get_3d_shifts([ad_am_1_05, ad_am_neg_1_04],
                                       pixel_scale=0.004,
                                       sub_pixel_fraction=sub_pix_frac)

        assert len(shifts["y_shifts"]) == approx(10 / sub_pix_frac, rel=0.2)


class TestGetImagingWaveset:
    def test_returns_default_wave_range_when_passed_no_effects(self):
        wave_bin_edges = fov_mgr.get_imaging_waveset([])
        assert len(wave_bin_edges) == 0

    def test_returns_waveset_of_filter(self):
        filt = _filter_tophat_curve()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fov_mgr.get_imaging_waveset([filt], **kwargs)
        assert len(wave_bin_edges) == 2

    def test_returns_waveset_of_psf(self):
        psf = _const_psf()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fov_mgr.get_imaging_waveset([psf], **kwargs)
        assert len(wave_bin_edges) == 3

    def test_returns_waveset_of_psf_and_filter(self):
        filt = _filter_tophat_curve()
        psf = _const_psf()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fov_mgr.get_imaging_waveset([filt, psf], **kwargs)
        assert len(wave_bin_edges) == 5

    def test_returns_waveset_of_ncpa_psf_inside_filter_edges(self):
        filt = _filter_tophat_curve()
        psf = _ncpa_psf()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fov_mgr.get_imaging_waveset([psf, filt], **kwargs)
        assert min(wave_bin_edges) == 1.
        assert max(wave_bin_edges) == 2.
        assert len(wave_bin_edges) == 9


class TestGetImagingHeaders:
    def test_returns_header_grid_detector_size(self):
        pass


class TestGetImagingFOVs:
    def test_returns_header_grid_detector_size(self):
        pass

