import pytest
from pytest import approx
from copy import deepcopy
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

from scopesim.optics import image_plane_utils as imp_utils
from scopesim.optics.fov_manager import chunk_edges
from scopesim.tests.mocks.py_objects import imagehdu_objects as imo


PLOTS = False


class TestSplitHeader:
    def test_is_the_header_split_into_the_right_amount_of_chunks(self):
        hdr = imp_utils.header_from_list_of_xy([-1.024, 1.024],
                                               [-1.024, 1.024], 0.002)
        hdrs = imp_utils.split_header(hdr, 128)
        area_sum = np.sum([hdr["NAXIS1"] * hdr["NAXIS2"] for hdr in hdrs])
        assert len(hdrs) == 64
        assert area_sum == hdr["NAXIS1"] * hdr["NAXIS2"]

    @pytest.mark.parametrize("x, y, pix", [(0.19, 0.2, 0.01),
                                           (2.19, 1.55, 0.01),
                                           (1.29, 1.2, 0.02),
                                           (2.55, 3.75, 1)])
    def test_final_header_is_smaller_for_odd_size(self, x, y, pix):
        hdr = imp_utils.header_from_list_of_xy([-1., x],
                                               [-1., y], pix)
        hdrs = imp_utils.split_header(hdr, 100)
        area_sum = np.sum([hdr["NAXIS1"] * hdr["NAXIS2"] for hdr in hdrs])
        assert area_sum == hdr["NAXIS1"] * hdr["NAXIS2"]


class TestAddImageHDUtoImageHDU:
    def big_small_hdus(self, big_wh=(20, 10), big_offsets=(0, 0),
                       small_wh=(6, 3), small_offsets=(0, 0), pixel_scale=0.1):
        w, h = np.array(big_wh) // 2
        x = np.array([-w, -w, w, w]) + big_offsets[0]
        y = np.array([h, -h, -h, h]) + big_offsets[1]
        big = imp_utils.header_from_list_of_xy(x, y, pixel_scale, "X")
        im = np.ones([big["NAXIS2"], big["NAXIS1"]])
        big = fits.ImageHDU(header=big, data=im)

        w, h = np.array(small_wh) // 2
        x = np.array([-w, -w, w, w]) + small_offsets[0]
        y = np.array([h, -h, -h, h]) + small_offsets[1]
        small = imp_utils.header_from_list_of_xy(x, y, pixel_scale, "X")
        im = np.ones([small["NAXIS2"], small["NAXIS1"]])
        small = fits.ImageHDU(header=small, data=im)

        return big, small

    def test_smaller_hdu_is_fully_in_larger_hdu(self):
        """yellow box in box"""
        big, small = self.big_small_hdus()
        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(small, big, wcs_suffix="X")

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == big_sum + small_sum

    def test_smaller_cube_is_fully_inside_larger_cube(self):
        """yellow box in box"""
        big, small = self.big_small_hdus()
        big.data = big.data[None, :, :] * np.ones(3)[:, None, None]
        small.data = small.data[None, :, :] * np.ones(3)[:, None, None]

        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(small, big, wcs_suffix="X")

        if PLOTS:
            plt.imshow(new.data[1, :, :], origin="lower")
            plt.show()

        assert np.sum(new.data) == big_sum + small_sum

    def test_larger_hdu_encompases_smaller_hdu(self):
        """monochrome box"""
        big, small = self.big_small_hdus()
        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(big, small, wcs_suffix="X")

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == 2 * small_sum

    def test_smaller_hdu_is_partially_in_larger_hdu(self):
        """yellow quarter top-right"""
        big, small = self.big_small_hdus(small_wh=(20, 10), small_offsets=(10, 5))
        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(small, big, wcs_suffix="X")

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == 1.25 * big_sum

    def test_larger_hdu_is_partially_in_smaller_hdu(self):
        """yellow quarter bottom-left"""
        big, small = self.big_small_hdus(small_wh=(20, 10), small_offsets=(10, 5))
        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(big, small, wcs_suffix="X")

        if PLOTS:

            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == 1.25 * big_sum

    def test_larger_cube_is_partially_in_smaller_cube(self):
        """yellow quarter bottom-left"""
        big, small = self.big_small_hdus(small_wh=(20, 10), small_offsets=(10, 5))
        big.data = big.data[None, :, :] * np.ones(3)[:, None, None]
        small.data = small.data[None, :, :] * np.ones(3)[:, None, None]

        big_sum, small_sum = np.sum(big.data), np.sum(small.data)
        new = imp_utils.add_imagehdu_to_imagehdu(big, small, wcs_suffix="X")

        if PLOTS:

            plt.imshow(new.data[1, :, :], origin="lower")
            plt.show()

        assert np.sum(new.data) == 1.25 * big_sum

    def test_larger_hdu_is_fully_outside_smaller_hdu(self):
        """monochrome box"""
        big, small = self.big_small_hdus(small_offsets=(15, 0))
        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(big, small, wcs_suffix="X")

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == small_sum


    # face is downloaded from the internet.
    @pytest.mark.webtest
    def test_python_image_coords(self):
        # numpy uses a system of im[y, x]
        # numpy.shape[0] = y_len, numpy.shape[1] = x_len
        # FITS uses bottom left as CRPIXn = [1, 1]
        # matplotlib just needs origin='lower' to display these correctly

        from scipy.datasets import face
        im = face()[::-1, :, 1]
        print(np.shape(face()))

        hdu = fits.ImageHDU(data=im)
        hdu.header["CDELT1"] = 1
        hdu.header["CDELT2"] = -1
        hdu.header["CRVAL1"] = 0
        hdu.header["CRVAL2"] = 0
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CTYPE1"] = "LINEAR"
        hdu.header["CTYPE2"] = "LINEAR"
        hdu.header["CRPIX1"] = im.shape[1]/2
        hdu.header["CRPIX2"] = im.shape[0]/2

        # hdu.writeto("racoon.fits", overwrite=True)
        if PLOTS:
            plt.imshow(hdu.data, origin="lower")
            plt.show()


class TestOverlayImage:
    def test_overlay_images_works_as_expected_for_2d_images(self):
        big = np.zeros((100, 100))
        small = np.ones((10, 10))
        im = imp_utils.overlay_image(small, big, (20, 80))
        assert np.sum(im[50:, :50]) == np.sum(small)

        if PLOTS:
            plt.imshow(im, origin="lower")
            plt.show()

    def test_overlay_images_works_for_3d_cubes(self):
        big = np.zeros((10, 100, 100))
        small = np.ones((10, 10, 10))
        im = imp_utils.overlay_image(small, big, (20, 80))
        assert np.sum(im) == np.sum(small)

        if PLOTS:
            plt.imshow(im[5, :, :], origin="lower")
            plt.show()

    def test_overlay_images_works_for_2d_images_at_border(self):
        big = np.zeros((100, 100))
        small = np.ones((10, 10))
        im = imp_utils.overlay_image(small, big, (0, 0))
        assert np.sum(im) == 0.25 * np.sum(small)

        if PLOTS:
            plt.imshow(im, origin="lower")
            plt.show()

    def test_overlay_images_works_for_3d_cubes_at_border(self):
        big = np.zeros((10, 100, 100))
        small = np.ones((10, 10, 10))
        im = imp_utils.overlay_image(small, big, (0, 0))
        assert np.sum(im) == 0.25 * np.sum(small)

        if PLOTS:
            plt.subplot(131)
            plt.imshow(im[:, :, 0], origin="lower")
            plt.subplot(132)
            plt.imshow(im[:, 0, :], origin="lower")
            plt.subplot(133)
            plt.imshow(im[0, :, :], origin="lower")
            plt.show()


class TestLatticeIsSharedBetweenChunks:
    """Regions carved from one origin in whole-pixel steps share a lattice.

    This is what chunked FOVs are: the image plane, the full-size chunks and
    the smaller last chunk are all measured from the same detector edge. When
    the region is not a whole number of pixels across, CRVAL used to absorb
    each region's own rounding remainder, putting them on lattices offset from
    one another by fractions of a pixel.
    """

    @pytest.mark.parametrize("extent_px", [24.0, 21.0, 17.32, 17.1, 12621 + 1/3])
    @pytest.mark.parametrize("chunk", [4, 8])
    def test_chunks_land_on_the_parent_lattice(self, extent_px, chunk):
        pixel_scale = 0.1
        lo = -extent_px / 2 * pixel_scale
        hi = lo + extent_px * pixel_scale

        parent, parent_naxis = imp_utils.create_wcs_from_points(
            np.array([[lo, lo], [hi, hi]]), pixel_scale)

        # the real chunk-edge helper, not a copy of it
        edges = [lo, *chunk_edges(lo, hi, chunk * pixel_scale), hi]

        for x0, x1 in zip(edges[:-1], edges[1:]):
            child, child_naxis = imp_utils.create_wcs_from_points(
                np.array([[x0, x0], [x1, x1]]), pixel_scale)
            # the child's pixel 0 must sit on a whole-pixel offset of the
            # parent's grid, otherwise the projection has to snap it
            world = child.wcs_pix2world([[0, 0]], 0)
            pix = parent.wcs_world2pix(world, 0)[0]
            np.testing.assert_allclose(pix, np.round(pix), atol=1e-6)

    @pytest.mark.parametrize("pnts, pixel_scale", [
        (np.array([[0, 0]]), 5),                                  # single point
        (np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]]), 5),      # sub-pixel
    ])
    def test_regions_smaller_than_a_pixel_stay_centred(self, pnts, pixel_scale):
        wcs, naxis = imp_utils.create_wcs_from_points(pnts, pixel_scale)
        np.testing.assert_array_equal(naxis, [1, 1])
        np.testing.assert_array_equal(wcs.wcs.crval,
                                      pnts.mean(axis=0))


class TestOverlayImageSnapping:
    """The origin, not the centre, must be snapped to the pixel grid.

    ``ceil(coords) - shape // 2`` snapped the centre instead, which is off by
    up to a whole pixel whenever ``coords`` is not on the lattice implied by
    ``small_im.shape``, in a direction that flips with the parity of that
    shape. Adjacent FOVs of different sizes were therefore displaced opposite
    ways, and their shared edge gained a duplicated or a dropped row/column.
    """

    @pytest.mark.parametrize("size", [1, 2, 3, 4, 7, 8])
    def test_places_image_on_its_own_lattice_without_shift(self, size):
        big = np.zeros((32, 32))
        small = np.ones((size, size))
        # centre coordinate that the shape's own lattice implies
        centre = 10 + (size - 1) / 2
        imp_utils.overlay_image(small, big, (centre, centre))
        rows = np.where(big.any(axis=1))[0]
        assert rows[0] == 10
        assert rows[-1] == 10 + size - 1

    @pytest.mark.parametrize("size", [1, 2, 3, 4, 7, 8])
    def test_snaps_sub_pixel_offset_the_same_way_for_any_shape(self, size):
        # A 0.3 pix offset must never move the image by a whole pixel, whatever
        # the parity of the shape.
        big = np.zeros((32, 32))
        small = np.ones((size, size))
        centre = 10 + (size - 1) / 2 + 0.3
        imp_utils.overlay_image(small, big, (centre, centre))
        rows = np.where(big.any(axis=1))[0]
        assert rows[0] == 10
        assert rows[-1] == 10 + size - 1

    def test_adjacent_differently_sized_images_tile_exactly(self):
        # Mimics chunked FOVs, where the last chunk is smaller than the rest.
        big = np.zeros((1, 20))
        origin = 0
        for size in (8, 8, 4):
            small = np.ones((1, size))
            imp_utils.overlay_image(small, big, (origin + (size - 1) / 2, 0))
            origin += size
        assert (big == 1).all()


class TestRescaleImageHDU:
    @pytest.mark.parametrize("pixel_scale", [0.3, 0.5, 1, 2, 3])
    def test_rescales_a_2D_imagehdu(self, pixel_scale):
        hdu0 = imo._image_hdu_rect()
        hdu1 = imp_utils.rescale_imagehdu(deepcopy(hdu0), pixel_scale)#/3600)

        hdr0 = hdu0.header
        hdr1 = hdu1.header

        assert hdr1["NAXIS1"] == np.ceil(hdr0["NAXIS1"] / pixel_scale)
        assert hdr1["NAXIS2"] == np.ceil(hdr0["NAXIS2"] / pixel_scale)

    @pytest.mark.parametrize("pixel_scale", [0.3, 0.5, 1, 2, 3])
    def test_rescales_a_3D_imagehdu(self, pixel_scale):
        hdu0 = imo._image_hdu_rect()
        hdu0.data = hdu0.data[None, :, :] * np.ones(5)[:, None, None]
        hdu1 = imp_utils.rescale_imagehdu(deepcopy(hdu0), pixel_scale)#/3600)

        hdr0 = hdu0.header
        hdr1 = hdu1.header

        assert np.sum(hdu0.data) == approx(np.sum(hdu1.data))
        assert hdr1["NAXIS1"] == np.ceil(hdr0["NAXIS1"] / pixel_scale)
        assert hdr1["NAXIS2"] == np.ceil(hdr0["NAXIS2"] / pixel_scale)
        assert hdr1["NAXIS3"] == hdr0["NAXIS3"]


class TestReorientImageHDU:
    @pytest.mark.parametrize("angle", [0, 30, -45])
    def test_reorients_a_2D_imagehdu(self, angle):
        hdu0 = imo._image_hdu_rect()
        angle *= np.pi / 180
        pc_matrix = {"PC1_1": np.cos(angle),
                     "PC1_2": np.sin(angle),
                     "PC2_1": -np.sin(angle),
                     "PC2_2": np.cos(angle)}

        hdu0.header.update(pc_matrix)
        hdu1 = imp_utils.reorient_imagehdu(deepcopy(hdu0))

        hdr0 = hdu0.header
        hdr1 = hdu1.header

        new_x = hdr0["NAXIS1"] * np.cos(abs(angle)) + \
                hdr0["NAXIS2"] * np.sin(abs(angle))
        new_y = hdr0["NAXIS1"] * np.sin(abs(angle)) + \
                hdr0["NAXIS2"] * np.cos(abs(angle))

        if PLOTS:
            plt.subplot(121)
            plt.imshow(hdu0.data, origin="lower")
            plt.subplot(122)
            plt.imshow(hdu1.data, origin="lower")
            plt.show()

        assert hdr1["NAXIS1"] == np.ceil(new_x)
        assert hdr1["NAXIS2"] == np.ceil(new_y)

    @pytest.mark.parametrize("angle", [0, 30, -45])
    def test_reorients_a_3D_imagehdu(self, angle):
        hdu0 = imo._image_hdu_rect()
        hdu0.data = hdu0.data[None, :, :] * np.ones(5)[:, None, None]

        angle *= np.pi / 180
        pc_matrix = {"PC1_1": np.cos(angle),
                     "PC1_2": np.sin(angle),
                     "PC2_1": -np.sin(angle),
                     "PC2_2": np.cos(angle)}

        hdu0.header.update(pc_matrix)
        hdu1 = imp_utils.reorient_imagehdu(deepcopy(hdu0))

        hdr0 = hdu0.header
        hdr1 = hdu1.header

        new_x = hdr0["NAXIS1"] * np.cos(abs(angle)) + \
                hdr0["NAXIS2"] * np.sin(abs(angle))
        new_y = hdr0["NAXIS1"] * np.sin(abs(angle)) + \
                hdr0["NAXIS2"] * np.cos(abs(angle))

        if PLOTS:
            plt.subplot(121)
            plt.imshow(hdu0.data[2,:,:], origin="lower")
            plt.subplot(122)
            plt.imshow(hdu1.data[2,:,:], origin="lower")
            plt.show()

        assert hdr1["NAXIS1"] == np.ceil(new_x)
        assert hdr1["NAXIS2"] == np.ceil(new_y)


class TestSubPixelFractions:
    @pytest.mark.parametrize("x, y, frac", [(0.5, 0.5, 0.25),
                                            (10.1, -9.1, 0.09)])
    def test_returns_expected_origin_fraction(self, x, y, frac):
        xs, ys, fracs = imp_utils.sub_pixel_fractions(x, y)
        assert fracs[0] == approx(frac)
        assert xs[0] == int(np.floor(x))
        assert ys[0] == int(np.floor(y))

    # def test_returns_expected_origin_fraction_for_list(self):
    #     x, y = np.array([1.1, 2.9]), np.array([0.0, 0.5])
    #     xs, ys, fracs = imp_utils.sub_pixel_fractions(x, y)
    #     print(xs)


class TestSkyDetWCS:
    def test_wcs_roundtrip(self):
        # FIXME: This should be a fixture, but everywhere in this module...
        hdu = imo._image_hdu_rect()
        sky_wcs = WCS(hdu.header)

        # Scale 1.0 from _image_hdu_rect cdelt, 20 arbitrary plate scale
        det_wcs, nax = imp_utils.det_wcs_from_sky_wcs(sky_wcs, 1.0, 20)
        new_wcs, nax = imp_utils.sky_wcs_from_det_wcs(det_wcs, 1.0, 20,
                                                      naxis=nax)

        # Compare header representations because of missing NAXIS info in new
        assert new_wcs.to_header() == sky_wcs.to_header()
        assert all(nax == (hdu.header["NAXIS1"], hdu.header["NAXIS2"]))


class TestWCSFromMinimalPoints:
    def test_all_zero_points(self):
        pnts = np.array([[0, 0]])
        wcs, naxis = imp_utils.create_wcs_from_points(pnts, 5)
        np.testing.assert_array_equal(naxis, [1, 1])
        np.testing.assert_array_equal(wcs.wcs.crpix, [1, 1])
        np.testing.assert_array_equal(wcs.wcs.crval, [0, 0])

    def test_all_within_one_pixel(self):
        pnts = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]])
        wcs, naxis = imp_utils.create_wcs_from_points(pnts, 5)
        np.testing.assert_array_equal(naxis, [1, 1])
        np.testing.assert_array_equal(wcs.wcs.crpix, [1, 1])
        np.testing.assert_array_equal(wcs.wcs.crval, [0, 0])

    def test_all_within_one_pixel_along_one_axis(self):
        pnts = np.array([[-3, 0], [-2, 0], [0, 0], [2, 0]])
        wcs, naxis = imp_utils.create_wcs_from_points(pnts, 1)
        np.testing.assert_array_equal(naxis, [5, 1])
        np.testing.assert_array_equal(wcs.wcs.crpix, [3, 1])
        np.testing.assert_array_equal(wcs.wcs.crval, [-0.5, 0])
