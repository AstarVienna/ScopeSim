import numpy as np
import pytest
from _pytest.python_api import approx
from matplotlib import pyplot as plt

from scopesim.base_classes import PoorMansHeader
from scopesim.effects import SpectralTraceList
from scopesim.optics import FieldOfView
from scopesim.optics import fov_manager_utils as fm_utils
from scopesim.optics.image_plane_utils import pix2val
from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim.tests.mocks.py_objects import aperture_objects as apo
from scopesim.tests.mocks.py_objects import trace_list_objects as tlo


@pytest.fixture(scope="function")
def full_trace_list():
    return tlo.make_trace_hdulist()


PLOTS = False


class TestGet3DShifts:
    def test_returns_zeros_when_no_shift3d_effects_passed(self):
        shifts = fm_utils.get_3d_shifts([], wave_min=0.7, wave_max=3,
                                        pixel_scale=0.004, sub_pixel_fraction=1)
        assert np.all(shifts["x_shifts"] == 0)
        assert np.all(shifts["y_shifts"] == 0)

    def test_returns_almost_zero_for_zenith_atmospheric_dispersion(self):
        ad_zenith = eo._atmospheric_dispersion(airmass=1.)
        shifts = fm_utils.get_3d_shifts([ad_zenith], pixel_scale=0.004,
                                        sub_pixel_fraction=1)
        assert np.all(shifts["x_shifts"] == 0)
        assert np.all(shifts["y_shifts"] == 0)

    def test_returns_non_zero_entries_for_off_zenith(self):
        ad_am_1_14 = eo._atmospheric_dispersion()
        shifts = fm_utils.get_3d_shifts([ad_am_1_14], pixel_scale=0.004,
                                        sub_pixel_fraction=1)
        assert np.all(shifts["x_shifts"] == 0)
        assert np.interp(0, shifts["y_shifts"][::-1],
                         shifts["wavelengths"][::-1]) == pytest.approx(1.5)

    def test_shift_transferred_to_x_axis_for_90_deg_pupil_angle(self):
        ad_am_1_05 = eo._atmospheric_dispersion(airmass=1.05, pupil_angle=90)
        shifts = fm_utils.get_3d_shifts([ad_am_1_05], pixel_scale=0.004,
                                        sub_pixel_fraction=1)
        assert np.interp(0, shifts["x_shifts"][::-1],
                         shifts["wavelengths"][::-1]) == pytest.approx(1.5)
        assert np.all(shifts["y_shifts"] == 0)

    def test_shifts_cancel_out_when_equal_and_opposite(self):
        ad_am_1_05 = eo._atmospheric_dispersion(airmass=1.05)
        ad_am_neg_1_05 = eo._atmospheric_dispersion(airmass=-1.05)
        shifts = fm_utils.get_3d_shifts([ad_am_1_05, ad_am_neg_1_05],
                                        pixel_scale=0.004, sub_pixel_fraction=1)

        assert len(shifts["y_shifts"]) == 2
        assert np.min(shifts["wavelengths"]) == ad_am_1_05.meta["wave_min"]
        assert np.max(shifts["wavelengths"]) == ad_am_1_05.meta["wave_max"]

    @pytest.mark.parametrize("sub_pix_frac", [1, 0.3, 0.1])
    def test_combined_shifts_reduced_to_usable_number(self, sub_pix_frac):
        ad_am_1_05 = eo._atmospheric_dispersion(airmass=1.05)
        ad_am_neg_1_04 = eo._atmospheric_dispersion(airmass=-1.04)
        shifts = fm_utils.get_3d_shifts([ad_am_1_05, ad_am_neg_1_04],
                                        pixel_scale=0.004,
                                        sub_pixel_fraction=sub_pix_frac)

        assert len(shifts["y_shifts"]) == approx(10 / sub_pix_frac, rel=0.2)


class TestGetImagingWaveset:
    def test_returns_default_wave_range_when_passed_no_effects(self):
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fm_utils.get_imaging_waveset([], **kwargs)
        assert len(wave_bin_edges) == 2

    def test_returns_waveset_of_filter(self):
        filt = eo._filter_tophat_curve()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fm_utils.get_imaging_waveset([filt], **kwargs)
        assert len(wave_bin_edges) == 2

    def test_returns_waveset_of_psf(self):
        psf = eo._const_psf()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fm_utils.get_imaging_waveset([psf], **kwargs)
        assert len(wave_bin_edges) == 4

    def test_returns_waveset_of_psf_and_filter(self):
        filt = eo._filter_tophat_curve()
        psf = eo._const_psf()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fm_utils.get_imaging_waveset([filt, psf], **kwargs)
        assert len(wave_bin_edges) == 4

    def test_returns_waveset_of_ncpa_psf_inside_filter_edges(self):
        filt = eo._filter_tophat_curve()
        psf = eo._ncpa_psf()
        kwargs = {"wave_min": 0.5, "wave_max": 2.5}
        wave_bin_edges = fm_utils.get_imaging_waveset([psf, filt], **kwargs)
        assert min(wave_bin_edges) == 1.
        assert max(wave_bin_edges) == 2.
        assert len(wave_bin_edges) == 9


class TestGetImagingHeaders:
    def test_throws_error_if_not_all_keywords_are_passed(self):
        apm = eo._img_aperture_mask()
        with pytest.raises(ValueError):
            fm_utils.get_imaging_headers([apm])

    def test_returns_set_of_headers_from_aperture_effects(self):
        apm = eo._img_aperture_mask(array_dict={"x": [-1.28, 1., 1., -1.28],
                                                "y": [-1.28, -1.28, 2., 2.]})
        kwargs = {"pixel_scale": 0.01, "plate_scale": 1,
                  "max_segment_size": 128**2, "chunk_size": 64}
        hdrs = fm_utils.get_imaging_headers([apm], **kwargs)

        area_sum = np.sum([hdr["NAXIS1"] * hdr["NAXIS2"] for hdr in hdrs])
        assert area_sum == 228 * 328

    def test_returns_set_of_headers_from_detector_list_effect(self):
        # det = eo._full_detector_list()
        det = eo._detector_list()
        kwargs = {"pixel_scale": 0.004, "plate_scale": 0.26666666666,
                  "max_segment_size": 2048 ** 2, "chunk_size": 1024}
        hdrs = fm_utils.get_imaging_headers([det], **kwargs)

        area_sum = np.sum([hdr["NAXIS1"] * hdr["NAXIS2"] for hdr in hdrs])
        assert area_sum == 4096**2

        if PLOTS:
            plt.subplot(121)
            for hdr in hdrs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(hdr)
                plt.plot(x*3600, y*3600)
                plt.title("Sky plane")
                plt.xlabel("[arcsec]")

            plt.subplot(122)
            for hdr in hdrs:
                from scopesim.optics.image_plane_utils import calc_footprint
                x, y = calc_footprint(hdr, "D")
                plt.plot(x, y)
                plt.title("Detector focal plane")
                plt.xlabel("[mm]")

            plt.show()


class TestGetImagingFOVs:
    def test_returns_fov_objects_for_basic_input(self):
        apm = eo._img_aperture_mask(array_dict={"x": [-1.0, 1.0, 1.0, -1.0],
                                                "y": [-1.0, -1.0, 1.0, 1.0]})
        kwargs = {"pixel_scale": 0.01, "plate_scale": 1,
                  "max_segment_size": 100 ** 2, "chunk_size": 100}

        hdrs = fm_utils.get_imaging_headers([apm], **kwargs)
        waveset = np.linspace(1, 2, 6)
        shifts = {"wavelengths": np.array([1, 2]),
                  "x_shifts": np.zeros(2),
                  "y_shifts": np.array([0, 1]) / 3600}  # 0..1 arcsec shift
        fovs = fm_utils.get_imaging_fovs(headers=hdrs, waveset=waveset,
                                         shifts=shifts)

        assert len(fovs) == (len(waveset)-1) * len(hdrs)

        if PLOTS:
            from scopesim.optics.image_plane_utils import calc_footprint
            plt.subplot(121)
            for fov in fovs:
                x, y = calc_footprint(fov.header)
                plt.fill(x*3600, y*3600, alpha=0.1, c="b")
                plt.title("Sky plane")
                plt.xlabel("[arcsec]")

            plt.subplot(122)
            for fov in fovs:
                x, y = calc_footprint(fov.header, "D")
                plt.fill(x, y)
                plt.title("Detector focal plane")
                plt.xlabel("[mm]")

            plt.show()


@pytest.mark.usefixtures("full_trace_list")
class TestGetSpectroscopyHeaders:
    def test_returns_headers(self, full_trace_list):
        params = {"pixel_scale": 0.1, "plate_scale": 0.1,
                  "wave_min": 0.7, "wave_max": 2.5}
        spt = SpectralTraceList(hdulist=full_trace_list, **params)
        apm = apo._basic_aperture()

        hdrs = fm_utils.get_spectroscopy_headers(effects=[spt, apm], **params)
        assert all([isinstance(hdr, PoorMansHeader) for hdr in hdrs])

        if PLOTS:
            for hdr in hdrs:
                x = np.array([0, hdr["NAXIS1"], hdr["NAXIS1"], 0])
                y = np.array([0, 0, hdr["NAXIS2"], hdr["NAXIS2"]])
                xw, yw = pix2val(hdr, x, y, "D")
                plt.plot(xw, yw)

            plt.show()


@pytest.mark.usefixtures("full_trace_list")
class TestGetSpectroscopyFOVs:
    def test_returns_fovs(self, full_trace_list):
        params = {"pixel_scale": 0.1, "plate_scale": 0.1,
                  "wave_min": 0.7, "wave_max": 2.5}
        spt = SpectralTraceList(hdulist=full_trace_list, **params)
        apm = apo._basic_aperture()

        shifts = {"wavelengths": np.array([0.7, 2.5]),
                  "x_shifts": np.array([0, 0]),
                  "y_shifts": np.array([0, 1/3600.])}
        hdrs = fm_utils.get_spectroscopy_headers(effects=[spt, apm], **params)
        fovs = fm_utils.get_spectroscopy_fovs(hdrs, shifts)

        assert all([isinstance(fov, FieldOfView) for fov in fovs])

        if PLOTS:
            ax1 = plt.subplot(121)
            ax2 = plt.subplot(122)
            for fov in fovs:
                sky, imp = fov.corners
                ax1.plot(sky[0], sky[1])
                ax2.plot(imp[0], imp[1])
            plt.show()
