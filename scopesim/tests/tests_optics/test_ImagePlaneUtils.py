import pytest
from pytest import approx
from copy import deepcopy
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from scopesim.optics import image_plane_utils as imp_utils
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

        # print([hdr["NAXIS1"] for hdr in hdrs], hdr["NAXIS1"])
        # print([hdr["NAXIS2"] for hdr in hdrs], hdr["NAXIS2"])


class TestAddImageHDUtoImageHDU:
    def big_small_hdus(self, big_wh=(20, 10), big_offsets=(0, 0),
                       small_wh=(6, 3), small_offsets=(0, 0), pixel_scale=0.1):
        w, h = np.array(big_wh) // 2
        x = np.array([-w, -w, w, w]) + big_offsets[0]
        y = np.array([h, -h, -h, h]) + big_offsets[1]
        big = imp_utils.header_from_list_of_xy(x, y, pixel_scale)
        im = np.ones([big["NAXIS2"], big["NAXIS1"]])
        big = fits.ImageHDU(header=big, data=im)

        w, h = np.array(small_wh) // 2
        x = np.array([-w, -w, w, w]) + small_offsets[0]
        y = np.array([h, -h, -h, h]) + small_offsets[1]
        small = imp_utils.header_from_list_of_xy(x, y, pixel_scale)
        im = np.ones([small["NAXIS2"], small["NAXIS1"]])
        small = fits.ImageHDU(header=small, data=im)

        return big, small

    def test_smaller_hdu_is_fully_in_larger_hdu(self):
        """yellow box in box"""
        big, small = self.big_small_hdus()
        big_sum, small_sum =  np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(small, big)

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == big_sum + small_sum

    def test_smaller_cube_is_fully_inside_larger_cube(self):
        """yellow box in box"""
        big, small = self.big_small_hdus()
        big.data = big.data[None, :, :] * np.ones(3)[:, None, None]
        small.data = small.data[None, :, :] * np.ones(3)[:, None, None]

        big_sum, small_sum =  np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(small, big)

        if PLOTS:
            plt.imshow(new.data[1, :, :], origin="lower")
            plt.show()

        assert np.sum(new.data) == big_sum + small_sum

    def test_larger_hdu_encompases_smaller_hdu(self):
        """monochrome box"""
        big, small = self.big_small_hdus()
        big_sum, small_sum =  np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(big, small)

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == 2 * small_sum

    def test_smaller_hdu_is_partially_in_larger_hdu(self):
        """yellow quarter top-right"""
        big, small = self.big_small_hdus(small_wh=(20, 10), small_offsets=(10, 5))
        big_sum, small_sum =  np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(small, big)

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == 1.25 * big_sum

    def test_larger_hdu_is_partially_in_smaller_hdu(self):
        """yellow quarter bottom-left"""
        big, small = self.big_small_hdus(small_wh=(20, 10), small_offsets=(10, 5))
        big_sum, small_sum =  np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(big, small)

        if PLOTS:

            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == 1.25 * big_sum

    def test_larger_cube_is_partially_in_smaller_cube(self):
        """yellow quarter bottom-left"""
        big, small = self.big_small_hdus(small_wh=(20, 10), small_offsets=(10, 5))
        big.data = big.data[None, :, :] * np.ones(3)[:, None, None]
        small.data = small.data[None, :, :] * np.ones(3)[:, None, None]

        big_sum, small_sum =  np.sum(big.data), np.sum(small.data)
        new = imp_utils.add_imagehdu_to_imagehdu(big, small)

        if PLOTS:

            plt.imshow(new.data[1, :, :], origin="lower")
            plt.show()

        assert np.sum(new.data) == 1.25 * big_sum

    def test_larger_hdu_is_fully_outside_smaller_hdu(self):
        """monochrome box"""
        big, small = self.big_small_hdus(small_offsets=(15, 0))
        big_sum, small_sum = np.sum(big.data), np.sum(small.data)

        new = imp_utils.add_imagehdu_to_imagehdu(big, small)

        if PLOTS:
            plt.imshow(new.data, origin="lower")
            plt.show()

        assert np.sum(new.data) == small_sum


    def test_python_image_coords(self):
        # numpy uses a system of im[y, x]
        # numpy.shape[0] = y_len, numpy.shape[1] = x_len
        # FITS uses bottom left as CRPIXn = [1, 1]
        # matplotlib just needs origin='lower' to display these correctly

        from scipy.misc import face
        im = face()[::-1, :, 1]
        print(np.shape(face()))

        hdu = fits.ImageHDU(data=im)
        hdu.header["CDELT1"] = 1
        hdu.header["CDELT2"] = -1
        hdu.header["CRVAL1"] = 0
        hdu.header["CRVAL2"] = 0
        hdu.header["CUNIT1"] = "DEG"
        hdu.header["CUNIT2"] = "DEG"
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


class TestRescaleImageHDU:
    @pytest.mark.parametrize("scale_factor", [0.3, 0.5, 1, 2, 3])
    def test_rescales_a_2D_imagehdu(self, scale_factor):
        hdu0 = imo._image_hdu_rect()
        hdu1 = imp_utils.rescale_imagehdu(deepcopy(hdu0), scale_factor/3600)

        hdr0 = hdu0.header
        hdr1 = hdu1.header

        assert hdr1["NAXIS1"] == np.ceil(hdr0["NAXIS1"] / scale_factor)
        assert hdr1["NAXIS2"] == np.ceil(hdr0["NAXIS2"] / scale_factor)

    @pytest.mark.parametrize("scale_factor", [0.3, 0.5, 1, 2, 3])
    def test_rescales_a_3D_imagehdu(self, scale_factor):
        hdu0 = imo._image_hdu_rect()
        hdu0.data = hdu0.data[None, :, :] * np.ones(5)[:, None, None]
        hdu1 = imp_utils.rescale_imagehdu(deepcopy(hdu0), scale_factor/3600)

        hdr0 = hdu0.header
        hdr1 = hdu1.header

        assert np.sum(hdu0.data) == approx(np.sum(hdu1.data))
        assert hdr1["NAXIS1"] == np.ceil(hdr0["NAXIS1"] / scale_factor)
        assert hdr1["NAXIS2"] == np.ceil(hdr0["NAXIS2"] / scale_factor)
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


