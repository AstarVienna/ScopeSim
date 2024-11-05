import numpy as np
from astropy import wcs
from astropy.io import fits


def _image_hdu_square(wcs_suffix=""):
    width = 100
    the_wcs = wcs.WCS(naxis=2, key=wcs_suffix)
    if wcs_suffix == "":
        #the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        the_wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        the_wcs.wcs.cunit = ["arcsec", "arcsec"]
    elif wcs_suffix == "D":
        the_wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        the_wcs.wcs.cunit = ["mm", "mm"]
    the_wcs.wcs.cdelt = [1, 1]
    the_wcs.wcs.crval = [0, 0]
    the_wcs.wcs.crpix = [(width + 1) / 2, (width + 1) / 2]

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
        # the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        the_wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        the_wcs.wcs.cunit = ["arcsec", "arcsec"]
    elif wcs_suffix == "D":
        the_wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        the_wcs.wcs.cunit = ["mm", "mm"]
    the_wcs.wcs.cdelt = [1, 1]
    the_wcs.wcs.crval = [0, 0]
    the_wcs.wcs.crpix = [(width + 1) / 2, (height + 1) / 2]
    the_wcs.wcs.pc = [[ca, sa], [-sa, ca]]

    image = np.random.random(size=(height, width))
    hdu = fits.ImageHDU(data=image, header=the_wcs.to_header())

    return hdu


def _image_hdu_three_wcs():
    hdu = _image_hdu_square()

    wcs_0 = wcs.WCS(hdu.header)
    wcs_d = wcs.WCS(naxis=2, key="D")
    wcs_g = wcs.WCS(naxis=2, key="G")

    wcs_d.wcs.ctype = ["LINEAR", "LINEAR"]
    wcs_g.wcs.ctype = ["LINEAR", "LINEAR"]

    wcs_d.wcs.cunit = ["mm", "mm"]
    wcs_d.wcs.crpix = wcs_0.wcs.crpix
    wcs_d.wcs.cdelt = [1., 1.]

    wcs_g.wcs.cunit = ["kpc", "kpc"]
    wcs_g.wcs.crpix = wcs_0.wcs.crpix
    wcs_g.wcs.cdelt = [1., 1.]
    hdu.header.update(wcs_0.to_header())
    hdu.header.update(wcs_d.to_header())
    hdu.header.update(wcs_g.to_header())

    return hdu
