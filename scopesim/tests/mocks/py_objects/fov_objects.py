import numpy as np
from astropy import units as u

from scopesim.optics import image_plane_utils as imp_utils
from scopesim.optics.fov import FieldOfView


def _centre_fov(n=55, waverange=(1.0, 2.0)):
    xsky = np.array([-n, n]) * u.arcsec.to(u.deg)
    ysky = np.array([-n, n]) * u.arcsec.to(u.deg)
    sky_hdr = imp_utils.header_from_list_of_xy(xsky, ysky, 1/3600.)
    imp_hdr = imp_utils.header_from_list_of_xy([-n, n], [-n, n], 1, "D")
    imp_hdr.update(sky_hdr)
    fov = FieldOfView(imp_hdr, waverange=waverange*u.um, area=1*u.m**2)

    return fov


def _centre_micado_fov(n=128, waverange=(1.9, 2.4)):
    """ n [arcsec] """
    xsky = np.array([-n, n]) * u.arcsec.to(u.deg)
    ysky = np.array([-n, n]) * u.arcsec.to(u.deg)
    pixscale = 0.004/3600.
    sky_hdr = imp_utils.header_from_list_of_xy(xsky, ysky, pixscale)
    imp_hdr = imp_utils.header_from_list_of_xy([-n, n], [-n, n], pixscale, "D")
    imp_hdr.update(sky_hdr)
    fov = FieldOfView(imp_hdr, waverange=waverange*u.um, area=1*u.m**2)

    return fov
