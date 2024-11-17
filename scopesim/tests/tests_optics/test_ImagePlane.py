"""Tests for ImagePlane and some ImagePlaneUtils"""

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from copy import deepcopy

import pytest
from pytest import approx

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from astropy import wcs

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import scopesim.optics.image_plane as opt_imp
import scopesim.optics.image_plane_utils as imp_utils

from scopesim.tests.mocks.py_objects.imagehdu_objects import \
    _image_hdu_square, _image_hdu_rect, _image_hdu_three_wcs,\
    _image_hdu_3d_data

PLOTS = False


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


@pytest.fixture(scope="function", name="input_table")
def fixture_input_table():
    x = [-10, -10, 0, 10, 10] * u.arcsec
    y = [-10, 10, 0, -10, 10] * u.arcsec
    f = [1, 3, 1, 1, 5]
    tbl = Table(names=["x", "y", "flux"], data=[x, y, f])

    return tbl


@pytest.fixture(scope="function", name="input_table_mm")
def fixture_input_table_mm():
    x = [-10, -10, 0, 10, 10] * u.mm
    y = [-10, 10, 0, -10, 10] * u.mm
    f = [1, 3, 1, 1, 5]
    tbl = Table(names=["x_mm", "y_mm", "flux"], data=[x, y, f])

    return tbl


class TestCombineTableBoundaries:
    def test_all_three_tables_are_inside_header_wcs(self, input_table):
        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)
        tbl3 = deepcopy(input_table)

        tbl2["x"] -= 25
        tbl3["y"] -= 60

        hdr = imp_utils._make_bounding_header_for_tables(tbl1, tbl2, tbl3)
        as2deg = u.arcsec.to(u.deg)
        for tbl in [tbl1, tbl2, tbl3]:
            x, y = imp_utils.val2pix(hdr, tbl["x"]*as2deg, tbl["y"]*as2deg)
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        if PLOTS:
            xy = imp_utils.calc_footprint(hdr)
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x, y)
            x0, y0 = imp_utils.val2pix(hdr, 0, 0)

            plt.plot(x, y, "b")
            plt.plot(x0, y0, "ro")
            for tbl in [tbl1, tbl2, tbl3]:
                x, y = imp_utils.val2pix(hdr, tbl["x"] / 3600.,
                                         tbl["y"] / 3600.)
                plt.plot(x, y, "k.")

            plt.show()

    def test_all_three_mm_tables_are_inside_header_wcs(self, input_table_mm):
        tbl1 = deepcopy(input_table_mm)
        tbl2 = deepcopy(input_table_mm)
        tbl3 = deepcopy(input_table_mm)

        tbl2["x_mm"] += 50
        tbl3["y_mm"] += 25

        hdr = imp_utils._make_bounding_header_for_tables(tbl1, tbl2, tbl3,
                                                         pixel_scale=100*u.um)
        for tbl in [tbl1, tbl2, tbl3]:
            x, y = imp_utils.val2pix(hdr, tbl["x_mm"], tbl["y_mm"], "D")
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        if PLOTS:
            xy = imp_utils.calc_footprint(hdr, "D")
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x, y, "D")
            x0, y0 = imp_utils.val2pix(hdr, 0, 0, "D")

            plt.plot(x, y, "b")
            plt.plot(x0, y0, "ro")
            for tbl in [tbl1, tbl2, tbl3]:
                x, y = imp_utils.val2pix(hdr, tbl["x_mm"], tbl["y_mm"], "D")
                plt.plot(x, y, "k.")

            plt.show()


class TestCombineImageHDUBoundaries:
    def test_all_two_imagehdus_are_inside_header_wcs(self, image_hdu_square,
                                                     image_hdu_rect):

        image_hdu_rect.header["CRVAL1"] -= 70  # * u.arcsec.to(u.deg)
        image_hdu_square.header["CRVAL2"] += 70  # * u.arcsec.to(u.deg)

        hdr = imp_utils._make_bounding_header_from_headers(
            image_hdu_square.header, image_hdu_rect.header)
        for imhdr in [image_hdu_square.header, image_hdu_rect.header]:
            xy = imp_utils.calc_footprint(imhdr)
            x, y = xy[:, 0], xy[:, 1]
            x *= u.arcsec.to(u.deg)
            y *= u.arcsec.to(u.deg)
            x, y = imp_utils.val2pix(hdr, x, y)
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        if PLOTS:
            for imhdr in [image_hdu_square.header, image_hdu_rect.header, hdr]:
                xy = imp_utils.calc_footprint(imhdr)
                x, y = xy[:, 0], xy[:, 1]
                xp, yp = imp_utils.val2pix(imhdr, x, y)
                plt.plot(x, y, "r-")

            plt.plot(0, 0, "ro")
            plt.gca().set_aspect(1)
            plt.show()

    def test_all_two_mm_imagehdus_are_inside_header_wcs(self,
                                                        image_hdu_square_mm,
                                                        image_hdu_rect_mm):

        image_hdu_rect_mm.header["CRVAL1D"] -= 40
        image_hdu_square_mm.header["CRVAL2D"] += 80

        hdr = imp_utils._make_bounding_header_from_headers(
            image_hdu_square_mm.header, image_hdu_rect_mm.header,
            pixel_scale=1*u.mm)
        for imhdr in [image_hdu_square_mm.header, image_hdu_rect_mm.header]:
            xy = imp_utils.calc_footprint(imhdr, "D")
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x, y, "D")
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        if PLOTS:
            for imhdr in [image_hdu_square_mm.header,
                          image_hdu_rect_mm.header, hdr]:
                xy = imp_utils.calc_footprint(imhdr, "D")
                x, y = xy[:, 0], xy[:, 1]
                xp, yp = imp_utils.val2pix(imhdr, x, y, "D")
                plt.plot(x, y, "r-")

            plt.plot(0, 0, "ro")
            plt.gca().set_aspect(1)
            plt.show()


class TestGetCanvasHeader:
    def test_all_5_objects_are_inside_header_wcs(self, image_hdu_square,
                                                 image_hdu_rect, input_table):

        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)
        tbl3 = deepcopy(input_table)

        tbl2["x"] -= 150
        tbl3["y"] -= 100

        image_hdu_rect.header["CRVAL1"] += 100 * u.arcsec.to(u.deg)
        image_hdu_square.header["CRVAL1"] += 0 * u.arcsec.to(u.deg)
        image_hdu_square.header["CRVAL2"] += 100 * u.arcsec.to(u.deg)

        hdr = imp_utils.get_canvas_header([image_hdu_square, tbl1, tbl2, tbl3,
                                           image_hdu_rect], pixel_scale=1*u.arcsec)

        as2deg = u.arcsec.to(u.deg)

        for im in [image_hdu_square.header, image_hdu_rect.header]:
            xy = imp_utils.calc_footprint(im)
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x*as2deg, y*as2deg)
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        for tbl in [tbl1, tbl2, tbl3]:
            x, y = imp_utils.val2pix(hdr, tbl["x"]*as2deg, tbl["y"]*as2deg)
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        if PLOTS:

            xy = imp_utils.calc_footprint(hdr)
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x, y)
            plt.plot(x, y, "b")
            x0, y0 = imp_utils.val2pix(hdr, 0, 0)
            plt.plot(x0, y0, "ro")

            for tbl in [tbl1, tbl2, tbl3]:
                x, y = imp_utils.val2pix(hdr, tbl["x"]*as2deg, tbl["y"]*as2deg)
                plt.plot(x, y, "k.")

            for im in [image_hdu_square.header, image_hdu_rect.header]:
                xy = imp_utils.calc_footprint(im)
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(hdr, x, y)
                plt.plot(x, y, "r-")

            plt.gca().set_aspect(1)
            plt.show()

    def test_all_5_objects_are_inside_mm_header_wcs(self, image_hdu_square_mm,
                                                    image_hdu_rect_mm,
                                                    input_table_mm):

        tbl1 = deepcopy(input_table_mm)
        tbl2 = deepcopy(input_table_mm)
        tbl3 = deepcopy(input_table_mm)

        tbl2["x_mm"] -= 150
        tbl3["y_mm"] -= 100

        image_hdu_rect_mm.header["CRVAL1D"] += 100
        image_hdu_square_mm.header["CRVAL1D"] += 0
        image_hdu_square_mm.header["CRVAL2D"] += 100

        hdr = imp_utils.get_canvas_header([image_hdu_square_mm, tbl1, tbl2,
                                           tbl3, image_hdu_rect_mm],
                                          pixel_scale=1*u.mm)

        for im in [image_hdu_square_mm.header, image_hdu_rect_mm.header]:
            xy = imp_utils.calc_footprint(im, "D")
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x, y, "D")
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        for tbl in [tbl1, tbl2, tbl3]:
            x, y = imp_utils.val2pix(hdr, tbl["x_mm"], tbl["y_mm"], "D")
            for xi, yi in zip(x, y):
                assert 0 <= xi < hdr["NAXIS1"]
                assert 0 <= yi < hdr["NAXIS2"]

        if PLOTS:

            xy = imp_utils.calc_footprint(hdr, "D")
            x, y = xy[:, 0], xy[:, 1]
            x, y = imp_utils.val2pix(hdr, x, y, "D")
            plt.plot(x, y, "b")
            x0, y0 = imp_utils.val2pix(hdr, 0, 0, "D")
            plt.plot(x0, y0, "ro")

            for tbl in [tbl1, tbl2, tbl3]:
                x, y = imp_utils.val2pix(hdr, tbl["x_mm"], tbl["y_mm"], "D")
                plt.plot(x, y, "k.")

            for im in [image_hdu_square_mm.header, image_hdu_rect_mm.header]:
                xy = imp_utils.calc_footprint(im, "D")
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(hdr, x, y, "D")
                plt.plot(x, y, "r-")

            plt.gca().set_aspect(1)
            plt.show()


class TestAddTableToImageHDU:
    def test_points_are_added_to_small_canvas(self, input_table):
        tbl1 = deepcopy(input_table)
        hdr = imp_utils.get_canvas_header([tbl1])

        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)
        canvas_hdu = imp_utils.add_table_to_imagehdu(tbl1, canvas_hdu)

        assert np.sum(canvas_hdu.data) == np.sum(tbl1["flux"])

        if PLOTS:
            # "top left is green, top right is yellow"
            plt.imshow(canvas_hdu.data, origin="lower")
            plt.show()

    def test_mm_points_are_added_to_small_canvas(self, input_table_mm):
        tbl1 = deepcopy(input_table_mm)
        hdr = imp_utils.get_canvas_header([tbl1], pixel_scale=2*u.mm)

        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)
        canvas_hdu = imp_utils.add_table_to_imagehdu(tbl1, canvas_hdu,
                                                     wcs_suffix="D")

        assert np.sum(canvas_hdu.data) == np.sum(tbl1["flux"])

        if PLOTS:
            # "top left is green, top right is yellow"
            plt.imshow(canvas_hdu.data, origin="lower")
            plt.show()

    # add a paramaterised test here to check border cases
    def test_points_are_added_to_massive_canvas(self, input_table):
        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)
        tbl3 = deepcopy(input_table)

        tbl1["y"] += 50
        tbl2["x"] += 20
        tbl3["x"] -= 25
        tbl3["y"] -= 25

        hdr = imp_utils.get_canvas_header([tbl1, tbl2, tbl3],
                                          pixel_scale=1*u.arcsec)

        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)

        for tbl in [tbl1, tbl2, tbl3]:
            canvas_hdu = imp_utils.add_table_to_imagehdu(tbl, canvas_hdu, True)

        total_flux = np.sum([tbl1["flux"], tbl2["flux"], tbl3["flux"]])
        assert np.sum(canvas_hdu.data) == total_flux

        if PLOTS:
            x, y = imp_utils.val2pix(hdr, 0, 0)
            plt.plot(x, y, "ro")
            "top left is green, top right is yellow"
            plt.imshow(canvas_hdu.data, origin="lower")
            plt.show()

    def test_mm_points_are_added_to_massive_canvas(self, input_table_mm):
        tbl1 = deepcopy(input_table_mm)
        tbl2 = deepcopy(input_table_mm)
        tbl3 = deepcopy(input_table_mm)

        tbl1["y_mm"] += 50
        tbl2["x_mm"] += 20
        tbl3["x_mm"] -= 25
        tbl3["y_mm"] -= 25

        hdr = imp_utils.get_canvas_header([tbl1, tbl2, tbl3],
                                          pixel_scale=1*u.mm)
        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)

        for tbl in [tbl1, tbl2, tbl3]:
            canvas_hdu = imp_utils.add_table_to_imagehdu(tbl, canvas_hdu,
                                                         True, wcs_suffix="D")

        total_flux = np.sum([tbl1["flux"], tbl1["flux"], tbl1["flux"]])
        assert np.sum(canvas_hdu.data) == total_flux

        if PLOTS:
            x, y = imp_utils.val2pix(hdr, 0, 0, "D")
            plt.plot(x, y, "ro")
            # "top left is green, top right is yellow"
            plt.imshow(canvas_hdu.data, origin="lower")
            plt.show()


class TestAddImageHDUToImageHDU:
    def test_image_is_added_to_small_canvas(self, image_hdu_rect,
                                            image_hdu_square):
        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1"] -= 150*u.arcsec.to(u.deg)
        im_hdu.header["CRVAL2"] += 40*u.arcsec.to(u.deg)
        hdr = imp_utils.get_canvas_header([im_hdu, image_hdu_square])

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

    def test_mm_image_is_added_to_small_canvas(self, image_hdu_rect_mm,
                                               image_hdu_square_mm):
        im_hdu = image_hdu_rect_mm
        im_hdu.header["CRVAL1D"] -= 150
        im_hdu.header["CRVAL2D"] += 40
        hdr = imp_utils.get_canvas_header([im_hdu, image_hdu_square_mm],
                                          pixel_scale=1*u.mm)

        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(im_hdu, canvas_hdu,
                                                        wcs_suffix="D")

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

    def test_image_and_tables_on_large_canvas(self, input_table, image_hdu_rect,
                                              image_hdu_square):
        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)

        tbl1["y"] -= 100
        tbl2["x"] += 100
        tbl2["y"] += 100

        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1"] -= 150  # *u.arcsec.to(u.deg)
        im_hdu.header["CRVAL2"] += 20  # *u.arcsec.to(u.deg)

        total_flux = (np.sum(tbl1["flux"]) + np.sum(tbl2["flux"]) +
                      im_hdu.data.sum() + image_hdu_square.data.sum())

        hdr = imp_utils.get_canvas_header([im_hdu, image_hdu_square,
                                           tbl1, tbl2], pixel_scale=3*u.arcsec)
        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)

        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(im_hdu, canvas_hdu)
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(image_hdu_square,
                                                        canvas_hdu)
        canvas_hdu = imp_utils.add_table_to_imagehdu(tbl1, canvas_hdu)
        canvas_hdu = imp_utils.add_table_to_imagehdu(tbl2, canvas_hdu)

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

    def test_mm_image_and_tables_on_large_canvas(self, input_table_mm,
                                                 image_hdu_rect_mm,
                                                 image_hdu_square_mm):
        image_hdu_rect = image_hdu_rect_mm
        image_hdu_square = image_hdu_square_mm
        tbl1 = deepcopy(input_table_mm)
        tbl2 = deepcopy(input_table_mm)

        tbl1["y_mm"] -= 100
        tbl2["x_mm"] += 100
        tbl2["y_mm"] += 100

        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1D"] -= 150
        im_hdu.header["CRVAL2D"] += 20

        hdr = imp_utils.get_canvas_header([im_hdu, image_hdu_square,
                                           tbl1, tbl2], pixel_scale=3*u.mm)
        im = np.zeros((hdr["NAXIS2"], hdr["NAXIS1"]))
        canvas_hdu = fits.ImageHDU(header=hdr, data=im)

        canvas_hdu = imp_utils.add_table_to_imagehdu(tbl1, canvas_hdu,
                                                     wcs_suffix="D")
        canvas_hdu = imp_utils.add_table_to_imagehdu(tbl2, canvas_hdu,
                                                     wcs_suffix="D")
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(im_hdu, canvas_hdu,
                                                        wcs_suffix="D")
        canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(image_hdu_square,
                                                        canvas_hdu,
                                                        wcs_suffix="D")

        total_flux = np.sum(tbl1["flux"]) + np.sum(tbl2["flux"]) + \
                     np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
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
    def test_add_many_tables_and_imagehdus(self, input_table, image_hdu_rect,
                                           image_hdu_square):
        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)

        tbl1["y"] -= 50
        tbl2["x"] += 50
        tbl2["y"] += 50

        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1"] -= 150*u.arcsec.to(u.deg)
        im_hdu.header["CRVAL2"] += 20*u.arcsec.to(u.deg)

        fields = [im_hdu, tbl1, tbl2, image_hdu_square]
        hdr = imp_utils.get_canvas_header(fields, pixel_scale=1 * u.arcsec)

        implane = opt_imp.ImagePlane(hdr)
        implane.add(fields)

        total_flux = np.sum(tbl1["flux"]) + np.sum(tbl2["flux"]) + \
                     np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
        assert np.sum(implane.data) == approx(total_flux)

        if PLOTS:
            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im.header)
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(implane.header, x, y)
                plt.plot(x, y, "r-")

            for tbl in [tbl1, tbl2]:
                hdr = imp_utils._make_bounding_header_for_tables(tbl)
                xy = imp_utils.calc_footprint(hdr)
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(implane.header, x, y)
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(implane.header, 0, 0)
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(implane.data, origin="lower", norm=LogNorm())
            plt.show()

    def test_add_many_mm_tables_and_imagehdus(self, input_table_mm,
                                              image_hdu_rect_mm,
                                              image_hdu_square_mm):
        image_hdu_rect = image_hdu_rect_mm
        image_hdu_square = image_hdu_square_mm
        tbl1 = deepcopy(input_table_mm)
        tbl2 = deepcopy(input_table_mm)

        tbl1["y_mm"] -= 50
        tbl2["x_mm"] += 50
        tbl2["y_mm"] += 50

        im_hdu = image_hdu_rect
        im_hdu.header["CRVAL1D"] -= 150 # mm
        im_hdu.header["CRVAL2D"] += 20

        fields = [im_hdu, tbl1, tbl2, image_hdu_square]
        hdr = imp_utils.get_canvas_header(fields, pixel_scale=1*u.mm)
        implane = opt_imp.ImagePlane(hdr)
        implane.add(fields, wcs_suffix="D")

        total_flux = np.sum(tbl1["flux"]) + np.sum(tbl2["flux"]) + \
                     np.sum(im_hdu.data) + np.sum(image_hdu_square.data)
        assert np.sum(implane.data) == approx(total_flux)

        if PLOTS:
            for im in [im_hdu, image_hdu_square]:
                xy = imp_utils.calc_footprint(im.header, "D")
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(implane.header, x, y, "D")
                plt.plot(x, y, "r-")

            for tbl in [tbl1, tbl2]:
                hdr = imp_utils._make_bounding_header_for_tables(tbl,
                                                                 pixel_scale=1*u.mm)
                xy = imp_utils.calc_footprint(hdr, "D")
                x, y = xy[:, 0], xy[:, 1]
                x, y = imp_utils.val2pix(implane.header, x, y, "D")
                plt.plot(x, y, "r-")

            x0, y0 = imp_utils.val2pix(implane.header, 0, 0, "D")
            plt.plot(x0, y0, "ro")
            plt.gca().set_aspect(1)

            plt.imshow(implane.data, origin="lower", norm=LogNorm())
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
# ..todo: When you have time, reintegrate these tests, There are some good ones

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
    def test_header_contains_future_naxis_pixel_sizes(self, image_hdu_square,
                                                      image_hdu_rect):
        hdr = imp_utils.get_canvas_header([image_hdu_square, image_hdu_rect])
        assert hdr["NAXIS1"] == 100 + 1
        assert hdr["NAXIS2"] == 200 + 1

    @pytest.mark.parametrize("offset", -np.random.randint(200, 1001, 10))
    def test_header_contains_spread_out_regions(self, offset, image_hdu_square,
                                                image_hdu_rect):
        image_hdu_rect.header["CRVAL1"] += offset  # *u.arcsec.to(u.deg)
        hdr = imp_utils.get_canvas_header([image_hdu_square, image_hdu_rect])
        image_width = image_hdu_square.header["NAXIS1"] // 2 + \
                      image_hdu_rect.header["NAXIS1"] // 2 + abs(offset) + 1

        assert hdr["NAXIS1"] == image_width

    @pytest.mark.parametrize("offset", [0, 10, 25, 50])
    def test_header_has_correct_size_based_on_table_extremes(self, offset,
                                                             input_table):
        tbl1 = input_table
        tbl2 = deepcopy(input_table)
        tbl3 = deepcopy(input_table)
        tbl2["x"] += offset
        tbl3["y"] += offset
        hdr = imp_utils.get_canvas_header([tbl1, tbl2, tbl3],
                                          pixel_scale=0.1*u.arcsec)

        assert hdr["NAXIS1"] == np.max(tbl1["x"] + tbl2["x"]) * 10 + 3  # 4
        assert hdr["NAXIS2"] == np.max(tbl1["y"] + tbl3["y"]) * 10 + 3  # 4

    # @pytest.mark.parametrize("pix_scl", [5, 1, 0.5, 0.1])
    # def test_header_has_correct_size_with_tbl_and_image_input(self, input_table,
    #                                                           image_hdu_square,
    #                                                           pix_scl):
    #     input_table["x"] += 100
    #     hdr = imp_utils.get_canvas_header([image_hdu_square, input_table],
    #                                       pixel_scale=pix_scl * u.arcsec)
    #     assert hdr["NAXIS1"] == approx(hdr["CRPIX1"] +
    #                                    np.max(input_table["x"]) / pix_scl + 1,
    #                                    abs=0.5)

    def test_header_does_not_contain_negative_naxis_keywords(self):
        pass

#
# class TestGetCornerSkyCoordsFromTable:
#     def test_table_with_column_units_returns_right_values(self):
#         x, y = [1] * u.arcmin, [1] * u.arcmin
#         tbl = Table(names=["x", "y"], data=[x, y])
#         xsky, ysky = imp_utils.get_corner_sky_coords_from_table(tbl)
#
#         assert xsky[0]*u.deg == x[0].to(u.deg)
#         assert ysky[0]*u.deg == y[0].to(u.deg)
#
#     def test_table_with_meta_units_returns_right_values(self):
#         x, y = [1], [1]
#         tbl = Table(names=["x", "y"], data=[x, y])
#         tbl.meta.update({"x_unit": u.arcmin, "y_unit": u.arcmin})
#         xsky, ysky = impl_utils.get_corner_sky_coords_from_table(tbl)
#
#         assert xsky[0] == x[0]*u.arcmin.to(u.deg)
#         assert ysky[0] == y[0]*u.arcmin.to(u.deg)
#
#     def test_table_with_default_units_returns_right_values(self):
#         x, y = [60], [60]      # because default unit is arcsec
#         tbl = Table(names=["x", "y"], data=[x, y])
#         xsky, ysky = impl_utils.get_corner_sky_coords_from_table(tbl)
#
#         assert pytest.approx(xsky[0] == x[0] * u.arcmin.to(u.deg))
#         assert pytest.approx(ysky[0] == y[0] * u.arcmin.to(u.deg))
#
#
# class TestGetCornerSkyCoords:
#     def test_returns_coords_for_combination_of_table_and_header(self,
#                                                             image_hdu_square):
#         x, y = [-100, 70] * u.arcsec, [0, 0] * u.arcsec
#         tbl = Table(names=["x", "y"], data=[x, y])
#         xsky, ysky = impl_utils.get_corner_sky_coords([tbl, image_hdu_square])
#
#         assert np.all((xsky == x.to(u.deg).value))
#


class TestAddTableToImageHDU:
    @pytest.mark.parametrize("xpix, ypix, value",
                             [(51, 51, 2),
                              (48, 51, 3),
                              (48, 48, 4),
                              (51, 48, 5)])
    def test_integer_pixel_fluxes_are_added_correctly(self, xpix, ypix, value,
                                                      image_hdu_square):
        # Given the weird behaviour on pixel boundaries
        x, y = [1.5, -1.5, -1.5, 1.5]*u.arcsec, [1.5, 1.5, -1.5, -1.5]*u.arcsec
        flux = [1, 2, 3, 4] * u.Unit("ph s-1")
        tbl = Table(names=["x", "y", "flux"], data=[x, y, flux])

        hdu = imp_utils.add_table_to_imagehdu(tbl, image_hdu_square,
                                              sub_pixel=False)

        assert hdu.data[ypix, xpix] == value

    @pytest.mark.parametrize("x, y, flux, xpix, ypix, value",
                             [([0], [0], [1], 50, 50, .25),
                              ([0.7], [0.7], [1], 50, 50, 0.64),
                              ([-0.7], [-0.7], [1], 48, 48, 0.04),
                              ([5.5], [-4.7], [1], 55, 45, 0.8),
                              ([5.5], [-4.7], [1], 55, 44, 0.2)])
    def test_sub_pixel_fluxes_are_added_correctly(self, x, y, flux, xpix, ypix,
                                                  value, image_hdu_square):
        # Given the weird behaviour on pixel boundaries
        tbl = Table(names=["x", "y", "flux"],
                    data=[x*u.arcsec, y*u.arcsec, flux*u.Unit("ph s-1")])
        hdu = imp_utils.add_table_to_imagehdu(tbl, image_hdu_square,
                                              sub_pixel=True)

        if PLOTS:
            plt.imshow(hdu.data[45:55, 45:55], origin="lower")
            plt.colorbar()
            plt.show()

        assert hdu.data[ypix, xpix] == approx(value + 1)

    @pytest.mark.parametrize("x, y, flux",
                             [([100, -100], [0, 0], [10, 10])])
    def test_source_outside_canvas_are_ignored(self, x, y, flux,
                                               image_hdu_square):
            # Given the weird behaviour on pixel boundaries
            tbl = Table(names=["x", "y", "flux"],
                        data=[x * u.arcsec, y * u.arcsec,
                              flux * u.Unit("ph s-1")])
            hdu = imp_utils.add_table_to_imagehdu(tbl, image_hdu_square,
                                                  sub_pixel=True)

            assert np.sum(hdu.data) == np.prod(hdu.data.shape)


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

    def test_initialises_with_header_with_hdu(self, image_hdu_square,
                                              image_hdu_rect):
        hdr = imp_utils.get_canvas_header(pixel_scale=0.1 * u.arcsec,
                                          hdu_or_table_list=[image_hdu_rect,
                                                             image_hdu_square])
        implane = opt_imp.ImagePlane(hdr)
        assert isinstance(implane, opt_imp.ImagePlane)
        assert isinstance(implane.hdu, fits.ImageHDU)

    def test_throws_error_if_header_does_not_have_valid_wcs(self):
        with pytest.raises(ValueError):
            opt_imp.ImagePlane(fits.Header())


class TestImagePlaneAdd:
    def test_simple_add_imagehdu_conserves_flux(self, image_hdu_square,
                                                image_hdu_rect):
        fields = [image_hdu_rect, image_hdu_square]
        hdr = imp_utils.get_canvas_header(pixel_scale=1 * u.arcsec,
                                          hdu_or_table_list=fields)

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

    def test_simple_add_table_conserves_flux(self, image_hdu_rect):
        x = [75, -75]*u.arcsec
        y = [0, 0]*u.arcsec
        flux = [30, 20] * u.Unit("ph s-1")
        tbl = Table(names=["x", "y", "flux"], data=[x, y, flux])

        hdr = imp_utils.get_canvas_header(pixel_scale=0.1 * u.arcsec,
                                          hdu_or_table_list=[image_hdu_rect,
                                                             tbl])
        implane = opt_imp.ImagePlane(hdr)
        implane.add(tbl)
        assert np.isclose(np.sum(implane.data), np.sum(flux.value))

    def test_compound_add_image_and_table_conserves_flux(self, image_hdu_rect):
        x = [75, -75]*u.arcsec
        y = [0, 0]*u.arcsec
        flux = [30, 20] * u.Unit("ph s-1")
        tbl = Table(names=["x", "y", "flux"], data=[x, y, flux])

        hdr = imp_utils.get_canvas_header(pixel_scale=0.1 * u.arcsec,
                                          hdu_or_table_list=[image_hdu_rect,
                                                             tbl])
        in_sum = np.sum(flux.value) + np.sum(image_hdu_rect.data)

        implane = opt_imp.ImagePlane(hdr)
        implane.add(tbl)
        implane.add(image_hdu_rect)
        out_sum = np.sum(implane.data)

        assert out_sum == approx(in_sum, rel=5e-3)
