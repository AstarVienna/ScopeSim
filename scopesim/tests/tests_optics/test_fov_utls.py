from pytest import approx
import pytest

from matplotlib import pyplot as plt
import numpy as np
from synphot import Empirical1D, SourceSpectrum
from synphot.units import PHOTLAM
from astropy import units as u

from scopesim.optics.fov import FieldOfView, extract_range_from_spectrum
from scopesim.optics import image_plane_utils as imp_utils

from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.tests.mocks.py_objects import source_objects as so

PLOTS = False


@pytest.fixture(scope="function")
def cube_source():
    return so._cube_source()


@pytest.fixture(scope="function")
def basic_fov_header():
    return ho._basic_fov_header()


class TestExtractAreaFromImageHDU:
    def test_returns_full_cube_for_thick_fov(self, cube_source,
                                             basic_fov_header):
        fov = FieldOfView(basic_fov_header, [0.5, 2.5])
        field = cube_source.fields[0]
        new_field = fov.extract_area_from_imagehdu(field, fov.get_corners("deg")[0])

        if PLOTS:
            xy = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(xy[:, 0], xy[:, 1], c="r", label="FOV")
            xy = imp_utils.calc_footprint(field.header)
            plt.fill(xy[:, 0], xy[:, 1], c="y", label="Source")
            xy = imp_utils.calc_footprint(new_field.header)
            plt.fill(xy[:, 0], xy[:, 1], c="g", label="Extracted")

            plt.grid()
            plt.legend()
            plt.gca().set_aspect("equal")
            plt.show()

        assert new_field.header["NAXIS1"] == field.header["NAXIS1"]
        assert new_field.header["NAXIS2"] == field.header["NAXIS2"]
        assert new_field.header["NAXIS3"] == field.header["NAXIS3"]

    def test_returns_wavelength_cut_cube_for_thin_fov(self, cube_source,
                                                      basic_fov_header):
        fov = FieldOfView(basic_fov_header, [1.3, 1.7])
        field = cube_source.fields[0]
        new_field = fov.extract_area_from_imagehdu(field, fov.get_corners("deg")[0])

        if PLOTS:
            xy = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(xy[:, 0], xy[:, 1], c="r", label="FOV")
            xy = imp_utils.calc_footprint(field.header)
            plt.fill(xy[:, 0], xy[:, 1], c="y", label="Source")
            xy = imp_utils.calc_footprint(new_field.header)
            plt.fill(xy[:, 0], xy[:, 1], c="g", label="Extracted")

            plt.grid()
            plt.legend()
            plt.gca().set_aspect("equal")
            plt.show()

        assert new_field.header["NAXIS1"] == field.header["NAXIS1"]
        assert new_field.header["NAXIS2"] == field.header["NAXIS2"]
        assert new_field.header["NAXIS3"] == 21

    def test_returns_eigth_cube_for_3d_offset_fov(self, cube_source,
                                                  basic_fov_header):
        hdr = basic_fov_header
        hdr["CRVAL1"] += 75 * hdr["CDELT1"]
        hdr["CRVAL2"] += 75 * hdr["CDELT2"]
        fov = FieldOfView(hdr, [1.5, 2.5])
        field = cube_source.fields[0]
        new_field = fov.extract_area_from_imagehdu(field, fov.get_corners("deg")[0])

        if PLOTS:
            xy = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(xy[:, 0], xy[:, 1], c="r", label="FOV")
            xy = imp_utils.calc_footprint(field.header)
            plt.fill(xy[:, 0], xy[:, 1], c="y", label="Source")
            xy = imp_utils.calc_footprint(new_field.header)
            plt.fill(xy[:, 0], xy[:, 1], c="g", label="Extracted")

            plt.grid()
            plt.legend()
            plt.gca().set_aspect("equal")
            plt.show()

        # Note: 26 is correct because there are actually 25.5 source pixels in
        #       the FOV, but the cutout is "generous".
        assert new_field.header["NAXIS1"] == 26
        assert new_field.header["NAXIS2"] == 26
        assert new_field.header["NAXIS3"] == 51


class TestExtractRangeFromSpectrum:
    def test_extracts_the_wave_range_needed(self):
        wave = np.arange(0.7, 2.5, 0.1) * u.um
        flux = np.arange(len(wave)) * PHOTLAM
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

        waverange = [1.98, 2.12] * u.um
        new_spec = extract_range_from_spectrum(spec, waverange)

        assert len(new_spec.waverange) == 2
        assert new_spec.waverange[0] == 1.98 * u.um
        assert new_spec(1.98 * u.um).value == approx(12.8)

    @pytest.mark.parametrize(("endpoint", "msg"),
                             [pytest.param(1.5, "Waverange does not overlap", marks=pytest.mark.xfail(reason="Check was disabled in function, dunno why.")),
                              (2.05, "Waverange only partially overlaps")])
    def test_logs_msg_for_waverang_overlap_mismatch(
            self, endpoint, msg, caplog):
        wave = np.arange(0.7, endpoint, 0.1) * u.um
        flux = np.arange(len(wave)) * PHOTLAM
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

        waverange = [1.98, 2.12] * u.um
        extract_range_from_spectrum(spec, waverange)

        assert msg in caplog.text
