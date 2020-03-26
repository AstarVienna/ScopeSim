import numpy as np
from astropy import wcs
from astropy.io import fits


def _image_hdu_square(wcs_suffix=""):
    width = 100
    the_wcs = wcs.WCS(naxis=2, key=wcs_suffix)
    if wcs_suffix == "":
        the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        the_wcs.wcs.cunit = ["arcsec", "arcsec"]
    elif wcs_suffix == "D":
        the_wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        the_wcs.wcs.cunit = ["mm", "mm"]
    the_wcs.wcs.cdelt = [1, 1]
    the_wcs.wcs.crval = [0, 0]
    the_wcs.wcs.crpix = [width // 2, width // 2]

    # theta = 24
    # ca, sa = np.cos(np.deg2rad(theta)), np.sin(np.deg2rad(theta))
    # the_wcs.wcs.pc = np.array([[ca, sa], [-sa, ca]])

    image = np.ones((width, width))
    hdu = fits.ImageHDU(data=image, header=the_wcs.to_header())

    return hdu


def _image_hdu_rect(wcs_suffix=""):
    width = 50
    height = 200
    angle = 0
    ca, sa = np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))
    the_wcs = wcs.WCS(naxis=2, key=wcs_suffix)
    if wcs_suffix == "":
        the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        the_wcs.wcs.cunit = ["arcsec", "arcsec"]
    elif wcs_suffix == "D":
        the_wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        the_wcs.wcs.cunit = ["mm", "mm"]
    the_wcs.wcs.cdelt = [1, 1]
    the_wcs.wcs.crval = [0, 0]
    the_wcs.wcs.crpix = [width // 2, height // 2]
    the_wcs.wcs.pc = [[ca, sa], [-sa, ca]]

    image = np.random.random(size=(height, width))
    hdu = fits.ImageHDU(data=image, header=the_wcs.to_header())

    return hdu
