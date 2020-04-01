import numpy as np
from astropy import wcs

from scopesim.optics import image_plane_utils as imp_utils
from scopesim.optics.image_plane_utils import header_from_list_of_xy


def _basic_fov_header():
    w, h = 150, 150
    skywcs = wcs.WCS(naxis=2)
    skywcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    skywcs.wcs.cdelt = [0.1, 0.1]
    skywcs.wcs.cunit = ["arcsec", "arcsec"]
    skywcs.wcs.crval = [0, 0]
    skywcs.wcs.crpix = [w / 2, h / 2]

    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [1, 1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    skyhdr = skywcs.to_header()
    dethdr = detwcs.to_header()
    skyhdr.update(dethdr)
    skyhdr["NAXIS"] = 2
    skyhdr["NAXIS1"] = w
    skyhdr["NAXIS2"] = h

    return skyhdr


def _implane_header():

    w, h = 150, 150
    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [1, 1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    dethdr = detwcs.to_header()
    dethdr["NAXIS"] = 2
    dethdr["NAXIS1"] = w
    dethdr["NAXIS2"] = h

    return dethdr


def _fov_header():
    w, h = 100, 100
    skywcs = wcs.WCS(naxis=2)
    skywcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    skywcs.wcs.cdelt = [0.2, 0.2]
    skywcs.wcs.cunit = ["arcsec", "arcsec"]
    skywcs.wcs.crval = [0, 0]
    skywcs.wcs.crpix = [w / 2, h / 2]

    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [1, 1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    skyhdr = skywcs.to_header()
    dethdr = detwcs.to_header()
    skyhdr.update(dethdr)
    skyhdr["NAXIS"] = 2
    skyhdr["NAXIS1"] = w
    skyhdr["NAXIS2"] = h

    return skyhdr


def _basic_dtcr_header(n=20, pix_size=0.01):
    xs = [-pix_size * n/2, pix_size * n/2]
    hdr = header_from_list_of_xy(xs, xs, pix_size, "D")
    return hdr


def _short_micado_slit_header():
    x = np.array([-1.5, 1.5]) / 3600.
    y = np.array([-0.01, 0.01]) / 3600.
    pix_scale_deg = 0.004 / 3600.
    header = imp_utils.header_from_list_of_xy(x, y, pix_scale_deg)
    header["APERTURE"] = 0

    return header


def _long_micado_slit_header():
    x = np.array([-1.5, 13.5]) / 3600.
    y = np.array([-0.01, 0.01]) / 3600.
    pix_scale_deg = 0.004 / 3600.
    header = imp_utils.header_from_list_of_xy(x, y, pix_scale_deg)
    header["APERTURE"] = 0

    return header








