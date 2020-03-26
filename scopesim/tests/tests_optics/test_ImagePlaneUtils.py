import pytest
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from scopesim.optics import image_plane_utils as imp_utils

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
    def test_add_one_image_to_another_in_the_right_position(self):

        big = imp_utils.header_from_list_of_xy(x=np.array([-20, -20, 20, 20]),
                                               y=np.array([10, -10, -10, 10]),
                                               pixel_scale=0.1)
        im = np.zeros([big["NAXIS2"], big["NAXIS1"]])
        big = fits.ImageHDU(header=big, data=im)

        small = imp_utils.header_from_list_of_xy(x=np.array([-3, -3, 3, 3])+10,
                                                 y=np.array([1, -1, -1, 1])-5,
                                                 pixel_scale=0.1)
        im = np.ones([small["NAXIS2"], small["NAXIS1"]])
        small = fits.ImageHDU(header=small, data=im)

        big = imp_utils.add_imagehdu_to_imagehdu(small, big)
        ycen, xcen = np.array(big.data.shape) // 2
        assert np.sum(big.data[:ycen, xcen:]) == np.sum(small.data)

        plt.imshow(big.data, origin="lower")
        plt.show()

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
    def test_overlay_images_works_as_expected(self):
        big = np.zeros((100, 100))
        small = np.ones((10, 10))
        im = imp_utils.overlay_image(small, big, (20, 80))
        assert np.sum(im[50:, :50]) == np.sum(small)

        if PLOTS:
            plt.imshow(im, origin="lower")
            plt.show()
