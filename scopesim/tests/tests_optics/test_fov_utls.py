from pytest import approx
import pytest

from matplotlib import pyplot as plt
import numpy as np
from synphot import Empirical1D, SourceSpectrum
from synphot.units import PHOTLAM
from astropy import units as u

from scopesim.optics import FieldOfView, fov_utils
from scopesim.optics import image_plane_utils as imp_utils

from scopesim.tests.mocks.py_objects.header_objects import _basic_fov_header
from scopesim.tests.mocks.py_objects.source_objects import _cube_source

PLOTS = False


@pytest.fixture(scope="function")
def cube_source():
    return _cube_source()


@pytest.fixture(scope="function")
def basic_fov_header():
    return _basic_fov_header()


@pytest.mark.usefixtures("cube_source", "basic_fov_header")
class TestExtractAreaFromImageHDU:
    def test_returns_full_cube_for_thick_fov(self, cube_source,
                                             basic_fov_header):
        fov = FieldOfView(basic_fov_header, [0.5, 2.5])
        field = cube_source.fields[0]
        new_field = fov_utils.extract_area_from_imagehdu(field, fov.volume())

        if PLOTS:
            x, y = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(x, y, c="r")
            x, y = imp_utils.calc_footprint(field.header)
            plt.fill(x, y, c="y")
            x, y = imp_utils.calc_footprint(new_field.header)
            plt.fill(x, y, c="g")

            plt.show()

        assert new_field.header["NAXIS1"] == field.header["NAXIS1"]
        assert new_field.header["NAXIS2"] == field.header["NAXIS2"]
        assert new_field.header["NAXIS3"] == field.header["NAXIS3"]

    def test_returns_wavelength_cut_cube_for_thin_fov(self, cube_source,
                                                      basic_fov_header):
        fov = FieldOfView(basic_fov_header, [1.3, 1.7])
        field = cube_source.fields[0]
        new_field = fov_utils.extract_area_from_imagehdu(field, fov.volume())

        if PLOTS:
            x, y = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(x, y, c="r")
            x, y = imp_utils.calc_footprint(field.header)
            plt.fill(x, y, c="y")
            x, y = imp_utils.calc_footprint(new_field.header)
            plt.fill(x, y, c="g")

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
        new_field = fov_utils.extract_area_from_imagehdu(field, fov.volume())

        if PLOTS:
            x, y = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(x, y, c="r")
            x, y = imp_utils.calc_footprint(field.header)
            plt.fill(x, y, c="y")
            x, y = imp_utils.calc_footprint(new_field.header)
            plt.fill(x, y, c="g")

            plt.show()

        assert new_field.header["NAXIS1"] == 26
        assert new_field.header["NAXIS2"] == 26
        assert new_field.header["NAXIS3"] == 51


class TestExtractRangeFromSpectrum:
    def test_extracts_the_wave_range_needed(self):
        wave = np.arange(0.7, 2.5, 0.1) * u.um
        flux = np.arange(len(wave)) * PHOTLAM
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

        waverange = [1.98, 2.12] * u.um
        new_spec = fov_utils.extract_range_from_spectrum(spec, waverange)

        assert len(new_spec.waverange) == 2
        assert new_spec.waverange[0] == 1.98 * u.um
        assert new_spec(1.98 * u.um).value == approx(12.8)

    def test_throws_error_if_no_overlap_between_waverange_and_waveset(self):
        wave = np.arange(0.7, 1.5, 0.1) * u.um
        flux = np.arange(len(wave)) * PHOTLAM
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

        with pytest.raises(ValueError):
            waverange = [1.98, 2.12] * u.um
            new_spec = fov_utils.extract_range_from_spectrum(spec, waverange)

    def test_throws_error_if_only_partial_overlap_exists(self):
        wave = np.arange(0.7, 2.05, 0.1) * u.um
        flux = np.arange(len(wave)) * PHOTLAM
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

        with pytest.raises(ValueError):
            waverange = [1.98, 2.12] * u.um
            new_spec = fov_utils.extract_range_from_spectrum(spec, waverange)
