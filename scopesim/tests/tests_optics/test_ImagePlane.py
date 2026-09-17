# -*- coding: utf-8 -*-
"""Tests for ImagePlane and some ImagePlaneUtils"""

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from copy import deepcopy
from itertools import product

import pytest
from pytest import approx

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import matplotlib.pyplot as plt

import scopesim.optics.image_plane as opt_imp
import scopesim.optics.image_plane_utils as imp_utils

from scopesim.tests.mocks.py_objects.imagehdu_objects import (
    _image_hdu_square,
    _image_hdu_rect,
    _image_hdu_three_wcs,
    _image_hdu_3d_data,
)

PLOTS = False


# The three functions below were oringinally stored in image_plane_utils,
# but ultimately only used here, so I moved them here as helpers.
# TODO: Check if the same functionality can be achieved in another way.
# ScopeSim internally uses something else (imageplane headers from detectors),
# so tests should ideally mock that instead...
def get_canvas_header(hdu_or_table_list, pixel_scale=1*u.arcsec):
    """
    Generate a fits.Header with a WCS that covers everything in the FOV.

    Parameters
    ----------
    hdu_or_table_list : list
        A list of Tables and/or ImageHDU py_objects

    pixel_scale : astropy.Quantity
        [arcsec] The pixel scale of the projection. Default in 1 arcsec

    Returns
    -------
    header : fits.Header
        A Header containing a WCS and NAXISn values to build an ImageHDU

    """
    def _get_headers(hdus_or_tables):
        for hdu_or_table in hdus_or_tables:
            if isinstance(hdu_or_table, fits.ImageHDU):
                yield hdu_or_table.header
            elif isinstance(hdu_or_table, fits.Header):
                yield hdu_or_table
            else:
                raise TypeError(
                    "hdu_or_table_list may only contain fits.ImageHDU "
                    f"or fits.Header, found {type(hdu_or_table)}.")

    headers = list(_get_headers(hdu_or_table_list))

    hdr = _make_bounding_header_from_headers(*headers, pixel_scale=pixel_scale)
    return hdr


def _make_bounding_header_from_headers(*headers, pixel_scale):
    """
    Return a Header with WCS and NAXISn keywords bounding all input ImageHDUs.

    Parameters
    ----------
    headers : list of fits.ImageHDU
    pixel_scale : u.Quantity
        [arcsec]

    Returns
    -------
    hdr : fits.Header

    """
    def _get_unit_from_headers(*headers, wcs_suffix: str) -> str:
        unit = headers[0][f"CUNIT1{wcs_suffix}"].lower()
        assert all(header[f"CUNIT{i}{wcs_suffix}"].lower() == unit
                   for header, i in product(headers, range(1, 3))), \
            [(i, header[f"CUNIT{i}{wcs_suffix}"])
             for header, i in product(headers, range(1, 3))]
        return unit

    wcs_suffix = "D" if pixel_scale.unit.physical_type == "length" else ""
    unit = u.Unit(_get_unit_from_headers(*headers, wcs_suffix=wcs_suffix))

    if unit.physical_type == "angle":
        unit = "deg"
        pixel_scale = pixel_scale.to_value(u.deg)
    else:
        pixel_scale = pixel_scale.to_value(unit)

    extents = [imp_utils.calc_footprint(header, wcs_suffix, unit)
               for header in headers]
    pnts = np.vstack(extents)

    hdr = imp_utils.header_from_list_of_xy(
        pnts[:, 0], pnts[:, 1], pixel_scale, wcs_suffix)
    hdr["NAXIS1"] += 1
    hdr["NAXIS2"] += 1
    hdr[f"CRVAL1{wcs_suffix}"] -= 0.5 * pixel_scale
    hdr[f"CRVAL2{wcs_suffix}"] -= 0.5 * pixel_scale

    return hdr


@pytest.fixture(scope="function", name="image_hdu_rect")
def fixture_image_hdu_rect():
    return _image_hdu_rect()


@pytest.fixture(scope="function", name="image_hdu_rect_mm")
def fixture_image_hdu_rect_mm():
    return _image_hdu_rect("D")


@pytest.fixture(scope="function", name="image_hdu_square")
def fixture_image_hdu_square():
    return _image_hdu_square()


@pytest.fixture(scope="function", name="image_hdu_square_mm")
def fixture_image_hdu_square_mm():
    return _image_hdu_square("D")


@pytest.fixture(scope="function", name="image_hdu_three_wcs")
def fixture_image_hdu_three_wcs():
    return _image_hdu_three_wcs()


@pytest.fixture(scope="function", name="image_hdu_3d_data")
def fixture_image_hdu_3d_data():
    return _image_hdu_3d_data()


class TestAddImageHDUToImageHDU:
    def test_image_is_added_to_small_canvas(
        self, image_hdu_rect, image_hdu_square,
    ):
        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1"] -= 150*u.arcsec.to(u.deg)
        im_hdu.header["CRVAL2"] += 40*u.arcsec.to(u.deg)
        hdr = get_canvas_header([im_hdu, image_hdu_square])

        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(im_hdu, canvas_hdu)
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(image_hdu_square, canvas_hdu)

        flux = np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
        assert np.sum(canvas_hdu.data) == approx(flux, rel=1e-2)

        if PLOTS:
            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im.header)
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(canvas_hdu.header, x, y)
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(canvas_hdu.header, 0, 0)
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(canvas_hdu.data, origin="lower")

            plt.show()

    def test_mm_image_is_added_to_small_canvas(
        self, image_hdu_rect_mm, image_hdu_square_mm,
    ):
        im_hdu = image_hdu_rect_mm
        im_hdu.header["CRVAL1D"] -= 150
        im_hdu.header["CRVAL2D"] += 40
        hdr = get_canvas_header([im_hdu, image_hdu_square_mm], 1*u.mm)

        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(
            im_hdu, canvas_hdu, wcs_suffix="D",
        )

        assert np.sum(canvas_hdu.data) == approx(np.sum(im_hdu.data))

        if PLOTS:
            for im in [im_hdu, image_hdu_square_mm]:
                xy = imp_utils.calc_footprint(im.header, "D")
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(canvas_hdu.header, x, y, "D")
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(canvas_hdu.header, 0, 0, "D")
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(canvas_hdu.data, origin="lower")

            plt.show()

    def test_images_on_large_canvas(self, image_hdu_rect, image_hdu_square):
        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1"] -= 150  # *u.arcsec.to(u.deg)
        im_hdu.header["CRVAL2"] += 20  # *u.arcsec.to(u.deg)

        total_flux = im_hdu.data.sum() + image_hdu_square.data.sum()

        hdr = get_canvas_header([im_hdu, image_hdu_square], 3*u.arcsec)
        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)

        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(im_hdu, canvas_hdu)
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(
            image_hdu_square, canvas_hdu,
        )

        assert np.sum(canvas_hdu.data) == approx(total_flux, rel=1e-2)

        if PLOTS:

            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im)
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(canvas_hdu, x, y)
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(canvas_hdu, 0, 0)
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(canvas_hdu.data, origin="lower")

            plt.show()

    def test_mm_images_on_large_canvas(
        self, image_hdu_rect_mm, image_hdu_square_mm,
    ):
        image_hdu_rect = image_hdu_rect_mm
        image_hdu_square = image_hdu_square_mm

        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1D"] -= 150
        im_hdu.header["CRVAL2D"] += 20

        hdr = get_canvas_header([im_hdu, image_hdu_square], 3*u.mm)
        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)

        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(
            im_hdu, canvas_hdu, wcs_suffix="D",
        )
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(
            image_hdu_square, canvas_hdu, wcs_suffix="D",
        )

        total_flux = np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
        assert np.sum(canvas_hdu.data) == approx(total_flux, rel=5e-3)

        if PLOTS:

            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im, "D")
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(canvas_hdu, x, y, "D")
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(canvas_hdu, 0, 0, "D")
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(canvas_hdu.data, origin="lower")

            plt.show()


class TestImagePlaneAdd:
    def test_simple_add_imagehdu_conserves_flux(
        self, image_hdu_square, image_hdu_rect,
    ):
        hdr = get_canvas_header([image_hdu_rect, image_hdu_square])

        orig_sum = image_hdu_rect.data.sum()

        print(wcs.WCS(image_hdu_rect))
        print(wcs.WCS(hdr))

        implane = opt_imp.ImagePlane(hdr)
        implane.add(image_hdu_rect)

        if PLOTS:
            plt.imshow(image_hdu_rect.data)
            x, y = wcs.WCS(image_hdu_rect).wcs_world2pix(0, 0, 1)
            print(x, y)
            plt.plot(x, y, "ro")
            plt.show()

            plt.imshow(implane.data)
            x, y = wcs.WCS(image_hdu_rect).wcs_world2pix(0, 0, 1)
            print(x, y)
            plt.plot(x, y, "ro")
            plt.show()

        assert np.sum(implane.data) == approx(orig_sum, rel=1e-2)

    def test_add_many_tables_and_imagehdus(self, image_hdu_rect, image_hdu_square):
        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1"] -= 150*u.arcsec.to(u.deg)
        im_hdu.header["CRVAL2"] += 20*u.arcsec.to(u.deg)

        hdr = get_canvas_header([im_hdu, image_hdu_square])

        implane = opt_imp.ImagePlane(hdr)
        implane.add([im_hdu, image_hdu_square])

        total_flux = np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
        assert np.sum(implane.data) == approx(total_flux)

        if PLOTS:
            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im.header)
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(implane.header, x, y)
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(implane.header, 0, 0)
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(implane.data, origin="lower", norm="log")
            plt.show()

    def test_add_many_mm_imagehdus(
            self, image_hdu_rect_mm, image_hdu_square_mm,
        ):
        image_hdu_rect = image_hdu_rect_mm
        image_hdu_square = image_hdu_square_mm

        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1D"] -= 150 # mm
        im_hdu.header["CRVAL2D"] += 20

        hdr = get_canvas_header([im_hdu, image_hdu_square], 1*u.mm)
        implane = opt_imp.ImagePlane(hdr)
        implane.add([im_hdu, image_hdu_square], wcs_suffix="D")

        total_flux = np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
        assert np.sum(implane.data) == approx(total_flux)

        if PLOTS:
            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im.header, "D")
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(implane.header, x, y, "D")
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(implane.header, 0, 0, "D")
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(implane.data, origin="lower", norm="log")
            plt.show()


class TestReorientImageHDU:
    def test_flux_remains_constant(self, image_hdu_rect):
        orig_sum = np.sum(image_hdu_rect.data)
        new_hdu = imp_utils.reorient_imagehdu(image_hdu_rect)
        new_sum = np.sum(new_hdu.data)

        assert new_sum == approx(orig_sum)

    def test_mm_flux_remains_constant(self, image_hdu_rect_mm):
        orig_sum = np.sum(image_hdu_rect_mm.data)
        new_hdu = imp_utils.reorient_imagehdu(image_hdu_rect_mm, wcs_suffix="D")
        new_sum = np.sum(new_hdu.data)

        assert new_sum == approx(orig_sum)


class TestRescaleImageHDU:
    @pytest.mark.parametrize("pixel_scale", [0.1, 0.237, 1, 2])
    def test_flux_remains_constant(self, image_hdu_rect, pixel_scale):
        orig_sum = np.sum(image_hdu_rect.data)
        new_hdu = imp_utils.rescale_imagehdu(image_hdu_rect,
                                             pixel_scale) #*u.arcsec.to(u.deg))
        new_sum = np.sum(new_hdu.data)

        assert new_sum == approx(orig_sum)

    @pytest.mark.parametrize("pixel_scale", [0.1, 0.237, 1, 2])
    def test_mm_flux_remains_constant(self, image_hdu_rect_mm, pixel_scale):
        orig_sum = np.sum(image_hdu_rect_mm.data)
        new_hdu = imp_utils.rescale_imagehdu(image_hdu_rect_mm, pixel_scale,
                                             wcs_suffix="D")
        new_sum = np.sum(new_hdu.data)

        assert new_sum == approx(orig_sum)

    @pytest.mark.parametrize("pixel_scale", [0.1, 0.237, 1, 2])
    def test_wcs_cdelt_scaled_correctly(self, image_hdu_three_wcs, pixel_scale):
        wcses = wcs.find_all_wcs(image_hdu_three_wcs.header)
        # this relies on find_all_wcs() sorting suffixes alphabetically
        fact = pixel_scale / wcses[0].wcs.cdelt[0]

        new_hdu = imp_utils.rescale_imagehdu(image_hdu_three_wcs, pixel_scale)
        new_wcses = wcs.find_all_wcs(new_hdu.header)
        assert new_wcses[0].wcs.cdelt[0] == pixel_scale
        assert new_wcses[0].wcs.cdelt[0] / fact == approx(wcses[0].wcs.cdelt[0])
        assert new_wcses[0].wcs.cdelt[1] / fact == approx(wcses[0].wcs.cdelt[1])
        assert new_wcses[1].wcs.cdelt[0] / fact == approx(wcses[1].wcs.cdelt[0])
        assert new_wcses[1].wcs.cdelt[1] / fact == approx(wcses[1].wcs.cdelt[1])
        assert new_wcses[2].wcs.cdelt[0] / fact == approx(wcses[2].wcs.cdelt[0])
        assert new_wcses[2].wcs.cdelt[1] / fact == approx(wcses[2].wcs.cdelt[1])

    def test_rescale_works_on_nondefault_wcs(self, image_hdu_three_wcs):
        pixel_scale = 2 * u.cm
        new_hdu = imp_utils.rescale_imagehdu(image_hdu_three_wcs,
                                             pixel_scale, "D")
        assert new_hdu.header['CDELT1D'] == 20


    def test_rescale_works_on_3d_imageplane(self, image_hdu_3d_data):
        pixel_scale = 0.274
        wcses = wcs.find_all_wcs(image_hdu_3d_data.header)
        fact = pixel_scale / wcses[0].wcs.cdelt[0]

        new_hdu = imp_utils.rescale_imagehdu(image_hdu_3d_data, pixel_scale)
        new_wcses = wcs.find_all_wcs(new_hdu.header)

        assert new_wcses[0].wcs.cdelt[0] == pixel_scale
        assert new_wcses[0].wcs.cdelt[2] == wcses[0].wcs.cdelt[2]
        assert new_wcses[1].wcs.cdelt[1] / fact == approx(wcses[1].wcs.cdelt[1])


###############################################################################
# TODO: When you have time, reintegrate these tests, There are some good ones

class TestGetSpatialExtentOfHeader:
    def test_returns_right_sky_coords_from_known_coords(self, image_hdu_square):
        xy = imp_utils.calc_footprint(image_hdu_square.header)
        xsky, ysky = xy[:, 0], xy[:, 1]
        xsky = np.array(xsky)
        xsky[xsky > 180 ] -= 360
        # xsky = np.array(xsky)*u.deg.to(u.arcsec)
        # ysky = np.array(ysky)*u.deg.to(u.arcsec)
        dx = max(xsky) - min(xsky)
        dy = max(ysky) - min(ysky)
        assert dx == approx(image_hdu_square.header["NAXIS1"])
        assert dy == approx(image_hdu_square.header["NAXIS2"])


class TestMakeImagePlaneHeader:
    def test_header_contains_future_naxis_pixel_sizes(
        self, image_hdu_square, image_hdu_rect,
    ):
        hdr = get_canvas_header([image_hdu_square, image_hdu_rect])
        assert hdr["NAXIS1"] == 100 + 1
        assert hdr["NAXIS2"] == 200 + 1

    @pytest.mark.parametrize("offset", -np.random.randint(200, 1001, 10))
    def test_header_contains_spread_out_regions(
        self, offset, image_hdu_square, image_hdu_rect,
    ):
        image_hdu_rect.header["CRVAL1"] += offset  # *u.arcsec.to(u.deg)
        hdr = get_canvas_header([image_hdu_square, image_hdu_rect])
        image_width = image_hdu_square.header["NAXIS1"] // 2 + \
                      image_hdu_rect.header["NAXIS1"] // 2 + abs(offset) + 1

        assert hdr["NAXIS1"] == image_width


class TestAddImagehduToImageHDU:
    @pytest.mark.parametrize("angle", [0, 30, 45, 89])
    def test_image_added_conserves_flux(self, angle, image_hdu_square):
        canvas = deepcopy(image_hdu_square)
        canvas.data = np.zeros((200, 200))
        canvas.header["CRPIX1"] *= 2
        canvas.header["CRPIX2"] *= 2

        angle = np.deg2rad(angle)
        image_hdu_square.data = np.ones((100, 100))
        image_hdu_square.header["PC1_1"] = np.cos(angle)
        image_hdu_square.header["PC1_2"] = np.sin(angle)
        image_hdu_square.header["PC2_1"] = -np.sin(angle)
        image_hdu_square.header["PC2_2"] = np.cos(angle)

        canvas = imp_utils.add_imagehdu_to_imagehdu(image_hdu_square, canvas)
        assert np.isclose(np.sum(canvas.data), np.sum(image_hdu_square.data))


class TestSubPixelFractions:
    @pytest.mark.parametrize("x, y, xx_exp ,yy_exp, ff_exp",
     [(   0,    0, [ 0, 0,  0, 0], [ 0,  0, 0, 0], [  1.,    0,    0,    0]),
      ( 0.2,  0.2, [ 0, 1,  0, 1], [ 0,  0, 1, 1], [0.64, 0.16, 0.16, 0.04]),
      (-0.2, -0.2, [-1, 0, -1, 0], [-1, -1, 0, 0], [0.04, 0.16, 0.16, 0.64]),
      ( 0.2, -0.2, [ 0, 1,  0, 1], [-1, -1, 0, 0], [0.16, 0.04, 0.64, 0.16])])
    def test_fractions_come_out_correctly_for_mixed_offsets(self, x, y, xx_exp,
                                                            yy_exp, ff_exp):
        xx, yy, ff = imp_utils.sub_pixel_fractions(x, y)
        for aa, bb in [[xx, xx_exp], [yy, yy_exp], [ff, ff_exp]]:
            assert all([a == approx(b) for a, b in zip(aa, bb)])


class TestImagePlaneInit:
    def test_throws_error_when_initialised_with_nothing(self):
        with pytest.raises(TypeError):
            opt_imp.ImagePlane()

    def test_initialises_with_header_with_hdu(
        self, image_hdu_square, image_hdu_rect,
    ):
        hdr = get_canvas_header([image_hdu_rect, image_hdu_square])
        implane = opt_imp.ImagePlane(hdr)
        assert isinstance(implane, opt_imp.ImagePlane)
        assert isinstance(implane.hdu, fits.ImageHDU)

    def test_throws_error_if_header_does_not_have_valid_wcs(self):
        with pytest.raises(ValueError):
            opt_imp.ImagePlane(fits.Header())
