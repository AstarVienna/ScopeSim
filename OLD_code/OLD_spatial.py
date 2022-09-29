"""
TODO: Insert docstring
"""

###############################################################################
# PlaneEffect
#
# DESCRIPTION
# The PlaneEffect functions are used to simulate effects that occur on all
# spectral layers equally, for example, sky rotation, telescope jitter,
# distortion, etc.
# To do this in the most general way possible, a PlaneEffect function contains
# 3 planes representing the deviation in position from an ideal optical train.
# The values in each of the 3 planes represent the distance the pixel should
# move in the x and y directions, and a weighting value.
#
# Several functions generate the various effects that occur. For example:
# - Rotation
# - Distortion
# - Translation
# - FlatField
# Some PlaneEffects only need to act on the positions of the incoming photons,
# e.g. ADC, while others are applicable to the whole array, e.g. Distortion,
# Flat field. As each PlaneEffect
# - CoordEffect
# - ArrayEffect
#
# As each PlaneEffect
#
# Classes:
#  CoordEffect
#  ArrayEffect
#
# Subclasses:
#  Rotation(ArrayEffect)
#  Distortion(ArrayEffect)
#  FlatField(ArrayEffect)
#  Translation(CoordEffect)

#
# Methods:
#
#
#

#from copy import deepcopy   ## Not used (OC)

import numpy as np
import scipy.ndimage as spi
from scipy.signal import fftconvolve

from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.io import fits

import scopesim.effects.effects_utils
import scopesim.effects.shifts
from scopesim import utils

__all__ = ["tracking", "derotator", "wind_jitter", "adc_shift",
           "make_distortion_maps", "get_distorion_offsets"]


def _gaussian_dist(x, mu, sig):
    p = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return p / np.sum(p)

def _linear_dist(x):
    p = np.array([1.]*len(x))
    return p / np.sum(p)

def _line_blur(arr, shift, kernel="gaussian", angle=0):
    """
    Introduce a linear blur due to tracking error.

    Parameters
    ----------
    - arr: [2D array] the image
    - shift: [pixel] how many pixels the image has moved

    Optional parameters
    -------------------
    - kernel: 'gaussian' - shift is the FWHM of the blur, approximating a random
                           walk in tracking error
              'linear' - shift is the length of the tracking blur with all
                         positions weighted equally, approximating no tracking
    - angle: [deg] the angle between image up and the zenith

    """

    # sample the shift at least every half pixel
    n = max(int(2 * shift) + 1, 3)
    if kernel == "gaussian":
        dr = np.linspace(-3 * shift, 3 * shift, 6 * n)
        weight = _gaussian_dist(dr, 0, shift)
    else:
        dr = np.linspace(0, shift, n)
        weight = _linear_dist(dr)

    dx = np.cos(np.deg2rad(angle)) * dr
    dy = np.sin(np.deg2rad(angle)) * dr

    tmp_arr = np.zeros(arr.shape)
    for x, y, w in zip(dx, dy, weight):  # TODO: x, y unused? (OC)
        tmp_arr += spi.shift(arr, (x, y), order=1) * w

    return tmp_arr

def _rotate_blur(arr, angle, kernel="gaussian"):
    """
    Introduce a rotational blur due to derotator error.

    Parameters
    ==========
    - arr: [2D array] the image
    - angle: [deg] the angle or rotation

    Optional parameters
    ===================
    - kernel: 'gaussian' - angle is the FWHM of the blur, approximating a random
                           walk in tracking error
              'linear' - shift is the length of the tracking blur with all
                         positions weighted equally, approximating no tracking
    """

    ang_at_cnr_pix = np.rad2deg(np.arctan2(1, np.sqrt(2) * arr.shape[0] // 2))
    n = max(3, int(angle / ang_at_cnr_pix) + 1)

    if kernel == "gaussian":
        d_ang = np.linspace(-3 * angle, 3 * angle, max(2, 6*n))
        weight = _gaussian_dist(d_ang, 0, angle)
    else:
        d_ang = np.linspace(0, angle, n)
        weight = _linear_dist(d_ang)

    tmp_arr = np.zeros(arr.shape)
    for ang, w in zip(d_ang, weight):
        tmp_arr += spi.rotate(arr, ang, order=1, reshape=False) * w

    return tmp_arr


def rotate_blur(image, angle):
    """
    Rotates and coadds an image over a given angle

    ..todo:
        Replace the function _rotate_blur in scopesim.spatial with this one

    """

    edge_pixel_angle = (np.arctan2(1, (psf_s.shape[0] // 2))*u.rad).to(u.deg).value

    angles = [angle]
    while angles[-1] > edge_pixel_angle and len(angles) < 25: angles += [angles[-1] / 2.]

    image_rot = np.copy(image)
    for ang in angles: image_rot += spi.rotate(image_rot, ang, reshape=False, order=1)

    image_rot /= image_rot.sum()

    return image_rot


def tracking(arr, cmds):
    """
    A method to simulate tracking errors
    ===== Currently a place holder with minimum functionality =========
    !! TODO, work out the shift during the DIT for the object RA, DEC etc !!
    """
    if cmds["SCOPE_DRIFT_DISTANCE"] > 0.:
        pix_res = cmds["SIM_PIXEL_SCALE"] / cmds["SIM_OVERSAMPLING"]
        kernel = cmds["SCOPE_DRIFT_PROFILE"]
        shift = cmds["SCOPE_DRIFT_DISTANCE"] / pix_res

        return _line_blur(arr, shift, kernel=kernel, angle=0)
    else:
        return arr


def derotator(arr, cmds):
    """
    A method to simulate field rotation in case the derotator is <100% effective
    ===== Currently a place holder with minimum functionality =========
    !! TODO, work out the rotation during the DIT for the object RA, DEC etc !!
    """
    if cmds["INST_DEROT_PERFORMANCE"] < 100.:
        eff = 1. - (cmds["INST_DEROT_PERFORMANCE"] / 100.)
        kernel = cmds["INST_DEROT_PROFILE"]
        angle = eff * cmds["OBS_EXPTIME"] * 15 / 3600.

        edge_smear = angle / cmds.pix_res
        if edge_smear > 50:
            print("The smear at the detector edge is large:", edge_smear, "[px]")

        return _rotate_blur(arr, angle, kernel=kernel)
    else:
        return arr


def wind_jitter(arr, cmds):
    """
    A method to simulate wind jitter
    ===== Currently a place holder with minimum functionality =========
    !! TODO, get the read spectrum for wind jitter !!
    !! Add in an angle parameter for the ellipse   !!
    """
    pix_res = cmds["SIM_PIXEL_SCALE"] / cmds["SIM_OVERSAMPLING"]
    fwhm = cmds["SCOPE_JITTER_FWHM"] / pix_res
    n = (fwhm / 2.35)
    kernel = Gaussian2DKernel(n, mode="oversample")

    return fftconvolve(arr, kernel, mode="same")
    #return convolve_fft(arr, kernel, allow_huge=True)


def adc_shift(cmds):
    """Generates a list of x and y shifts from a tests_commands object"""

    para_angle = cmds["OBS_PARALLACTIC_ANGLE"]
    effectiveness = cmds["INST_ADC_PERFORMANCE"] / 100.

    ## get the angle shift for each slice
    zenith_distance = utils.airmass2zendist(cmds["ATMO_AIRMASS"])
    angle_shift = [scopesim.effects.shifts.atmospheric_refraction(lam,
                                                                  zenith_distance,
                                                                  cmds["ATMO_TEMPERATURE"],
                                                                  cmds["ATMO_REL_HUMIDITY"],
                                                                  cmds["ATMO_PRESSURE"],
                                                                  cmds["SCOPE_LATITUDE"],
                                                                  cmds["SCOPE_ALTITUDE"])
                   for lam in cmds.lam_bin_centers]

    ## convert angle shift into number of pixels
    ## pixel shifts are defined with respect to last slice
    rel_shift = (angle_shift - angle_shift[-1])
    if np.max(np.abs(rel_shift)) > 1000:
        raise ValueError("Pixel shifts too great (>1000), check units")

    ## Rotate by the paralytic angle
    x = -rel_shift * np.sin(np.deg2rad(para_angle)) * (1. - effectiveness)
    y = -rel_shift * np.cos(np.deg2rad(para_angle)) * (1. - effectiveness)

    ## return values are in [arcsec]
    return x, y


def make_distortion_maps(real_xy, detector_xy, step=1):
    """
    Generate distortion maps based on star positions.

    The centres of the returned images correspond to the centre of the detector plane

    Parameters
    ----------
    real_xy : list
        [arcsec] Contains 2 arrays: ([x_pos], [y_pos]) where x_pos, y_pos are the
        coordinates of the real position of the stars

    detector_xy : list
        [arcsec] Contains 2 arrays: ([x_pos], [y_pos]) where x_pos, y_pos are the
        coordinates of the detected position of the stars

    step : float
        [arcsec] the grid spacing of the returned images

    Returns
    -------
    dx, dy : 2D array
        Returns two arrays with
    """

    from scipy.interpolate import griddata
    from astropy.io import fits

    dx = real_xy[0] - detector_xy[0]
    dy = real_xy[1] - detector_xy[1]

    w, h = np.max(np.abs(detector_xy[0])), np.max(np.abs(detector_xy[1]))
    gx, gy = np.arange(-w, w+1E-3, step),  np.arange(-h, h+1E-3, step)

    xx, yy = np.meshgrid(gx, gy)

    zx = griddata((real_xy[0], real_xy[1]), dx, (xx.flatten(), yy.flatten()),
                  method='cubic', fill_value=0).reshape(xx.shape)
    zy = griddata((real_xy[0], real_xy[1]), dy, (xx.flatten(), yy.flatten()),
                  method='cubic', fill_value=0).reshape(xx.shape)

    # the dx extension
    dx_hdu = fits.ImageHDU(data=zx.astype(np.float32))
    dx_hdu.header["CRVAL1"] = (0, "x dist from focal plane centre")
    dx_hdu.header["CRVAL2"] = (0, "y dist from focal plane centre")
    dx_hdu.header["CRPIX1"] = (w, "x coord ref pixel")
    dx_hdu.header["CRPIX2"] = (h, "y coord ref pixel")
    dx_hdu.header["CDELT1"] = (step, "x coord grid spacing")
    dx_hdu.header["CDELT2"] = (step, "y coord grid spacing")

    # the dy extension
    dy_hdu = fits.ImageHDU(data=zy.astype(np.float32),
                           header=dx_hdu.header)

    dx_hdu.header["COMMENT"] = "The distortion in the X direction"
    dy_hdu.header["COMMENT"] = "The distortion in the Y direction"


    tb_hdu = fits.BinTableHDU.from_columns([fits.Column(name='id',    format='E', array=np.arange(len(real_xy[0]))),
                                           fits.Column(name='x_real', format='E', array=real_xy[0]),
                                           fits.Column(name='y_real', format='E', array=real_xy[1]),
                                           fits.Column(name='x_obs',  format='E', array=real_xy[0]),
                                           fits.Column(name='y_obs',  format='E', array=real_xy[1]),
                                           fits.Column(name='del_x',  format='E', array=dx),
                                           fits.Column(name='del_y',  format='E', array=dy)])
    tb_hdu.header["COMMENT"] = "All units are [arcsec] relative to the FoV centre"

    # the primary header
    pri_hdu = fits.PrimaryHDU(header=dx_hdu.header)
    pri_hdu.header["EXT1"] = ("DX", "The x coordinate corrections")
    pri_hdu.header["EXT2"] = ("DY", "The y coordinate corrections")
    pri_hdu.header["EXT3"] = ("BinTable", "The data used to generate the distortion maps")

    # put it all together
    hdulist = fits.HDUList([pri_hdu, dx_hdu, dy_hdu, tb_hdu])

    return hdulist


def get_distorion_offsets(x, y, dist_map_hdus, corners):
    """
    Returns the distortion offsets for position relative to the FoV centre



    Parameters
    ----------
    x, y : float, list
        [arcsec] Distances from the centre of the distortion map (which should
        also be the centre of the FoV, given by CPREFn header Keywords)

    dist_map_hdus : list, astropy.io.fits.HDUList
        The distortion maps in X and Y directions. Accepts either
        - a list of two PrimaryHDU py_objects, each containing a map of the
          distortion on the x and y direction.
          E.g. ``dist_map_hdus=(hdu_dx, hdu_dy)`` where ``hdu_dx`` and ``hdu_dy``
          are FITS image py_objects (``ImageHDU``), or
        - a HDULlist object which contains 3 HDU py_objects: A PrimaryHDU and two
          ImageHDUs.
          - extension [0] (the PrimaryHDU) is nothing by a header,
          - extension [1] contains a map of the x-axis distortion,
          - extension [2] contains a map of the y-axis distortions

    corners : list
        [arcsec] A list containing 4 values for the borders of the distortion
        map grid in the following order (x_min, x_max, y_min, y_max)


    Returns
    -------
    dx, dy : float, list
        [arcsec] the shifts which need to be applied to the input positions x,y

    See Also
    --------
    make_distortion_maps()
    astropy.io.fits.HDUList, astropy.io.fits.PrimaryHDU

    """

    import scipy.interpolate as spi

    if isinstance(dist_map_hdus, fits.HDUList):
        xdist_hdu = dist_map_hdus[1]
        ydist_hdu = dist_map_hdus[2]
    elif isinstance(dist_map_hdus, (list, tuple)):
        xdist_hdu = dist_map_hdus[0]
        ydist_hdu = dist_map_hdus[1]


    if xdist_hdu.shape != ydist_hdu.shape:
        raise ValueError("Shape of X and Y distortion maps must be equal: "+str(xdist.shape)+" "+str(ydist.shape))

    xbins = np.linspace(corners[0], corners[1], xdist_hdu.header["NAXIS1"])
    ybins = np.linspace(corners[2], corners[3], xdist_hdu.header["NAXIS2"])

    xmesh, ymesh = np.meshgrid(xbins, ybins)
    xmesh, ymesh = xmesh.flatten(), ymesh.flatten()
    zmesh_dx = xdist_hdu.data.flatten()
    zmesh_dy = ydist_hdu.data.flatten()

    dx  = spi.griddata((xmesh, ymesh), zmesh_dx, (x, y), method="cubic", fill_value=0)
    dy = spi.griddata((xmesh, ymesh), zmesh_dy, (x, y), method="cubic", fill_value=0)

    return dx, dy
