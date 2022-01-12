import numpy as np
from scipy import ndimage as spi
from scipy.interpolate import RectBivariateSpline, griddata
from scipy.ndimage import zoom
from astropy import units as u
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from .. import rc, utils
from ..optics import image_plane_utils as imp_utils


def round_kernel_edges(kernel):
    y, x = np.array(kernel.shape).astype(int) // 2
    threshold = np.min([kernel[y, 0], kernel[y, -1],
                        kernel[0, x], kernel[-1, x]])
    kernel[kernel < threshold] = 0.

    return kernel


def nmrms_from_strehl_and_wavelength(strehl, wavelength, strehl_hdu,
                                     plot=False):
    """
    Returns the wavefront error needed to make a PSF with a desired strehl ratio

    Parameters
    ----------
    strehl : float
        [0.001, 1] Desired strehl ratio. Values 1<sr<100 will be scale to <1
    wavelength : float
        [um]
    strehl_hdu : np.ndarray
        2D map of strehl ratio as a function of wavelength [um] and residual
        wavefront error [nm RMS]
    plot : bool

    Returns
    -------
    nm : float
        [nm] residual wavefront error for generating an on-axis AnisoCADO PSF
        with the desired strehl ratio at a given wavelength

    """

    if 1. < strehl < 100.:
        strehl *= 0.01

    nm0 = strehl_hdu.header["CRVAL1"]
    dnm = strehl_hdu.header["CDELT1"]
    nm1 = strehl_hdu.header["NAXIS1"] * dnm + nm0
    nms = np.arange(nm0, nm1, dnm)

    w0 = strehl_hdu.header["CRVAL2"]
    dw = strehl_hdu.header["CDELT2"]
    w1 = strehl_hdu.header["NAXIS2"] * dw + w0
    ws = np.arange(w0, w1, dw)

    strehl_map = strehl_hdu.data

    nms_spline = RectBivariateSpline(ws, nms, strehl_map, kx=1, ky=1)
    strehls = nms_spline(wavelength, nms)[0]

    if strehl > np.max(strehls):
        raise ValueError("Strehl ratio ({}) is impossible at this wavelength "
                         "({}). Maximum Strehl possible is {}."
                         "".format(strehl, wavelength, np.max(strehls)))

    if strehls[0] < strehls[-1]:
        nm = np.interp(strehl, strehls, nms)
    else:
        nm = np.interp(strehl, strehls[::-1], nms[::-1])

    if plot:
        plt.plot(nms, strehls)
        plt.plot(nm, strehl, "ro")
        plt.show()

    return nm


def make_strehl_map_from_table(tbl, pixel_scale=1*u.arcsec):


    # pixel_scale = utils.quantify(pixel_scale, u.um).to(u.deg)
    # coords = np.array([tbl["x"], tbl["y"]]).T
    #
    # xmin, xmax = np.min(tbl["x"]), np.max(tbl["x"])
    # ymin, ymax = np.min(tbl["y"]), np.max(tbl["y"])
    # mesh = np.array(np.meshgrid(np.arange(xmin, xmax, pixel_scale),
    #                             np.arange(np.min(tbl["y"]), np.max(tbl["y"]))))
    # map = griddata(coords, tbl["layer"], mesh, method="nearest")
    #

    map = griddata(np.array([tbl.data["x"], tbl.data["y"]]).T,
                   tbl.data["layer"],
                   np.array(np.meshgrid(np.arange(-25, 26),
                                        np.arange(-25, 26))).T,
                   method="nearest")

    hdr = imp_utils.header_from_list_of_xy(np.array([-25, 25]) / 3600.,
                                           np.array([-25, 25]) / 3600.,
                                           pixel_scale=1/3600)

    map_hdu = fits.ImageHDU(header=hdr, data=map)

    return map_hdu


def rescale_kernel(image, scale_factor, spline_order=None):
    if spline_order is None:
        spline_order = utils.from_currsys("!SIM.computing.spline_order")
    sum_image = np.sum(image)
    image = zoom(image, scale_factor, order=spline_order)
    image = np.nan_to_num(image, copy=False)        # numpy version >=1.13

    # Re-centre kernel
    im_shape = image.shape
    dy, dx = np.divmod(np.argmax(image), im_shape[1]) - np.array(im_shape) // 2
    if dy > 0:
        image = image[2*dy:, :]
    elif dy < 0:
        image = image[:2*abs(dy), :]
    if dx > 0:
        image = image[:, 2*dx:]
    elif dx < 0:
        image = image[:, 2*abs(dx):]

    sum_new_image = np.sum(image)
    image *= sum_image / sum_new_image

    return image


def cutout_kernel(image, fov_header):
    h, w = image.shape
    xcen, ycen = 0.5 * w, 0.5 * h
    dx = 0.5 * fov_header["NAXIS1"]
    dy = 0.5 * fov_header["NAXIS2"]
    x0, x1 = max(0, int(xcen-dx)), min(w, int(xcen+dx))
    y0, y1 = max(0, int(ycen-dy)), min(w, int(ycen+dy))
    image_cutout = image[y0:y1, x0:x1]

    return image_cutout


def get_strehl_cutout(fov_header, strehl_imagehdu):

    image = np.zeros((fov_header["NAXIS2"], fov_header["NAXIS1"]))
    canvas_hdu = fits.ImageHDU(header=fov_header, data=image)
    canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(strehl_imagehdu,
                                                    canvas_hdu, spline_order=0,
                                                    conserve_flux=False)
    canvas_hdu.data = canvas_hdu.data.astype(int)

    return canvas_hdu


def nearest_index(x, x_array):
    # return int(round(np.interp(x, x_array, np.arange(len(x_array)))))
    return np.argmin(abs(x_array - x))


def get_psf_wave_exts(hdu_list, wave_key="WAVE0"):
    """
    Returns a dict of {extension : wavelength}

    Parameters
    ----------
    hdu_list

    Returns
    -------
    wave_set, wave_ext

    """

    if not isinstance(hdu_list, fits.HDUList):
        raise ValueError("psf_effect must be a PSF object: {}"
                         "".format(type(hdu_list)))

    tmp = np.array([[ii, hdu.header[wave_key]]
                    for ii, hdu in enumerate(hdu_list)
                    if wave_key in hdu.header and hdu.data is not None])
    wave_ext = tmp[:, 0].astype(int)
    wave_set = tmp[:, 1]

    # ..todo:: implement a way of getting the units from WAVEUNIT
    # until then assume everything is in um
    wave_set = utils.quantify(wave_set, u.um)

    return wave_set, wave_ext


def get_total_wfe_from_table(tbl):
    wfes = utils.quantity_from_table("wfe_rms", tbl, "um")
    n_surfs = tbl["n_surfaces"]
    total_wfe = np.sum(n_surfs * wfes**2)**0.5

    return total_wfe


def wfe2gauss(wfe, wave, width=None):
    strehl = wfe2strehl(wfe, wave)
    sigma = strehl2sigma(strehl)
    if width is None:
        width = int(np.ceil(8 * sigma))
        width += (width + 1) % 2
    gauss = sigma2gauss(sigma, x_size=width, y_size=width)

    return gauss


def wfe2strehl(wfe, wave):
    wave = utils.quantify(wave, u.um)
    wfe = utils.quantify(wfe, u.um)
    x = 2 * 3.1415926526 * wfe / wave
    strehl = np.exp(-x**2)
    return strehl


def strehl2sigma(strehl):
    amplitudes = [0.00465, 0.00480, 0.00506, 0.00553, 0.00637, 0.00793, 0.01092,
                  0.01669, 0.02736, 0.04584, 0.07656, 0.12639, 0.20474, 0.32156,
                  0.48097, 0.66895, 0.84376, 0.95514, 0.99437, 0.99982, 0.99999]
    sigmas = [19.9526, 15.3108, 11.7489, 9.01571, 6.91830, 5.30884, 4.07380,
              3.12607, 2.39883, 1.84077, 1.41253, 1.08392, 0.83176, 0.63826,
              0.48977, 0.37583, 0.28840, 0.22130, 0.16982, 0.13031, 0.1]
    sigma = np.interp(strehl, amplitudes, sigmas)
    return sigma


def sigma2gauss(sigma, x_size=15, y_size=15):
    kernel = Gaussian2DKernel(sigma, x_size=x_size, y_size=y_size,
                              mode="oversample").array
    kernel /= np.sum(kernel)
    return kernel


def rotational_blur(image, angle):
    """
    Rotates and coadds an image over a given angle to imitate a blur

    Parameters
    ----------
    image : array
        Image to blur
    angle : float
        [deg] Angle over which the image should be rotationally blurred

    Returns
    -------
    image_rot : array
        Blurred image

    """
    image_rot = np.copy(image)

    n_angles = 0
    rad_to_deg = 57.29578
    edge_pixel_unit_angle = np.arctan2(1, (image.shape[0] // 2)) * rad_to_deg
    while abs(angle) > edge_pixel_unit_angle and n_angles < 25:
        angle /= 2.
        image_rot += spi.rotate(image_rot, angle, reshape=False, order=1)
        n_angles += 1

    return image_rot / n_angles

def get_bkg_level(obj, bg_w):
    """
    Determine the background level of image or cube slices

    Returns a scalar if obj is a 2d image or a vector if obj is a 3D cube (one
    value for each plane).
    The method for background determination is decided by self.meta["bkg_width"]:
    If 0, the background is returned as zero (implying no background subtraction).
    If -1, the background is estimated as the median of the entire image (or
    cube plane).
    If positive, the background is estimated as the median of a frame of width
    `bkg_width` around the edges.
    """

    if obj.ndim == 2:
        if bg_w == 0:
            bkg_level = 0
        else:
            mask = np.zeros_like(obj, dtype=np.bool8)
            if bg_w > 0:
                mask[bg_w:-bg_w,bg_w:-bg_w] = True
            bkg_level = np.ma.median(np.ma.masked_array(obj, mask=mask))

    elif obj.ndim == 3:
        if bg_w == 0:
            bkg_level = np.array([0] * obj.shape[0])
        else:
            mask = np.zeros_like(obj, dtype=np.bool8)
            if bg_w > 0:
                mask[:, bg_w:-bg_w, bg_w:-bg_w] = True
            bkg_level = np.ma.median(np.ma.masked_array(obj, mask=mask),
                                     axis=(2, 1)).data

    else:
        raise ValueError("Unsupported dimension:", obj.ndim)
    return bkg_level
