from pytest import approx
import pytest

from matplotlib import pyplot as plt
import numpy as np
from synphot import Empirical1D, SourceSpectrum
from synphot.units import PHOTLAM
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

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

class TestExtractAreaKeepsTheInputPixelGrid:
    """The cutout is made of input pixels, so its WCS must say so.

    The cutout WCS used to be re-derived from the world corners *before* the
    slice was grown outwards to whole input pixels. That gave a WCS for a
    differently sized image (``NAXISn`` disagreed with ``data.shape``, silently
    corrected by ``fits.ImageHDU``, which left ``CRPIXn`` describing the wrong
    array) sitting on the FOV's pixel lattice rather than the input's. The
    cutout was then mis-registered by up to a whole pixel when projected, and
    by differing amounts for neighbouring FOVs, so their shared edge gained a
    duplicated or a dropped row/column.
    """

    @staticmethod
    def _image_hdu(naxis=32, pixel_scale=0.1, crval=(0.0, 0.0)):
        """Image with a distinct value in every pixel, in deg."""
        skywcs = WCS(naxis=2)
        skywcs.wcs.ctype = ["LINEAR", "LINEAR"]
        skywcs.wcs.cunit = ["deg", "deg"]
        skywcs.wcs.cdelt = 2 * [pixel_scale / 3600]
        skywcs.wcs.crval = list(crval)
        skywcs.wcs.crpix = 2 * [(naxis + 1) / 2]
        hdr = skywcs.to_header()
        hdr["BUNIT"] = ""
        data = np.arange(naxis * naxis, dtype=float).reshape(naxis, naxis)
        return fits.ImageHDU(header=hdr, data=data)

    @staticmethod
    def _fov(x_min, x_max, y_min, y_max, pixel_scale=0.1):
        skyhdr = imp_utils.header_from_list_of_xy(
            [x_min / 3600, x_max / 3600], [y_min / 3600, y_max / 3600],
            pixel_scale=pixel_scale / 3600)
        dethdr, _ = imp_utils.det_wcs_from_sky_wcs(
            WCS(skyhdr), pixel_scale, pixel_scale / 0.01)
        skyhdr.update(dethdr.to_header())
        return FieldOfView(skyhdr, [1.0, 2.0])

    # FOVs deliberately offset from the image's pixel edges by a fraction of a
    # pixel, which is what chunking produces whenever the detector extent is
    # not a whole number of pixels.
    @pytest.mark.parametrize("offset", [0.0, 0.16, 0.34, 0.5, 0.68, 0.84])
    def test_naxis_matches_data_shape(self, offset):
        hdu = self._image_hdu()
        fov = self._fov(-0.8 + offset * 0.1, 0.0 + offset * 0.1,
                        -0.8 + offset * 0.1, 0.0 + offset * 0.1)
        cut = fov.extract_area_from_imagehdu(hdu, fov.get_corners("deg")[0])
        assert (cut.header["NAXIS2"], cut.header["NAXIS1"]) == cut.data.shape

    @pytest.mark.parametrize("offset", [0.0, 0.16, 0.34, 0.5, 0.68, 0.84])
    def test_cutout_pixels_keep_their_world_coordinates(self, offset):
        hdu = self._image_hdu()
        fov = self._fov(-0.8 + offset * 0.1, 0.0 + offset * 0.1,
                        -0.8 + offset * 0.1, 0.0 + offset * 0.1)
        cut = fov.extract_area_from_imagehdu(hdu, fov.get_corners("deg")[0])

        # every cutout value must sit at the world coordinate its value had in
        # the input image
        in_wcs = WCS(hdu.header, naxis=2)
        cut_wcs = WCS(cut.header, naxis=2)
        ny, nx = cut.data.shape
        yy, xx = np.indices((ny, nx))
        world = cut_wcs.wcs_pix2world(
            np.column_stack([xx.ravel(), yy.ravel()]), 0)
        back = in_wcs.wcs_world2pix(world, 0).round(6)
        expected = hdu.data[back[:, 1].astype(int), back[:, 0].astype(int)]
        np.testing.assert_array_equal(cut.data.ravel(), expected)
        # ... and land on whole input pixels, not between them
        np.testing.assert_allclose(back, back.round(), atol=1e-6)
