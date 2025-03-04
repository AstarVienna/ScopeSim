# -*- coding: utf-8 -*-
"""Contains the base class for all PSF effects."""

from typing import ClassVar

import numpy as np
from scipy.signal import convolve
from scipy.ndimage import rotate
from astropy import units as u

from ..effects import Effect
from ...base_classes import ImagePlaneBase, FieldOfViewBase, FOVSetupBase
from ...utils import from_currsys, quantify, figure_factory, get_logger

logger = get_logger(__name__)

class PoorMansFOV:
    def __init__(self, pixel_scale, spec_dict, recursion_call=False):
        self.header = {"CDELT1": pixel_scale / 3600.,
                       "CDELT2": pixel_scale / 3600.,
                       "NAXIS": 2,
                       "NAXIS1": 128,
                       "NAXIS2": 128,
                       }
        self.meta = spec_dict
        self.wavelength = spec_dict["wave_mid"] * u.um
        if not recursion_call:
            self.hdu = PoorMansFOV(pixel_scale, recursion_call=True)


class PSF(Effect):
    z_order: ClassVar[tuple[int, ...]] = (40, 640)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        self.kernel = None
        self.valid_waverange = None
        self._waveset = []
        super().__init__(**kwargs)

        params = {
            "flux_accuracy": "!SIM.computing.flux_accuracy",
            "sub_pixel_flag": "!SIM.sub_pixel.flag",
            "convolve_mode": "same",      # "full", "same"
            "bkg_width": -1,
            "wave_key": "WAVE0",
            "normalise_kernel": True,
            "rounded_edges": True,
            "rotational_blur_angle": 0,
        }
        self.meta.update(params)
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta, self.cmds)
        self.convolution_classes = (FieldOfViewBase, ImagePlaneBase)

    def apply_to(self, obj, **kwargs):
        """Apply the PSF."""
        # 1. During setup of the FieldOfViews
        if isinstance(obj, FOVSetupBase) and self._waveset is not None:
            logger.debug("Executing %s, FoV setup", self.meta['name'])
            waveset = self._waveset
            if len(waveset) != 0:
                waveset_edges = 0.5 * (waveset[:-1] + waveset[1:])
                obj.split("wave", quantify(waveset_edges, u.um).value)

        # 2. During observe: convolution
        elif isinstance(obj, self.convolution_classes):
            logger.debug("Executing %s, convolution", self.meta['name'])
            if ((hasattr(obj, "fields") and len(obj.fields) > 0) or
                    (obj.hdu is not None)):
                kernel = self.get_kernel(obj).astype(float)

                # apply rotational blur for field-tracking observations
                rot_blur_angle = self.meta["rotational_blur_angle"]
                if abs(rot_blur_angle) > 0:
                    # makes a copy of kernel
                    kernel = rotational_blur(kernel, rot_blur_angle)

                # Round the edges of kernels so that the silly square stars
                # don't appear anymore
                if self.meta.get("rounded_edges", False) and kernel.ndim == 2:
                    kernel = self._round_kernel_edges(kernel)

                # normalise psf kernel      KERNEL SHOULD BE normalised within get_kernel()
                # if from_currsys(self.meta["normalise_kernel"], self.cmds):
                #    kernel /= np.sum(kernel)
                #    kernel[kernel < 0.] = 0.

                image = obj.hdu.data.astype(float)

                # subtract background level before convolving, re-add afterwards
                bkg_level = get_bkg_level(image, self.meta["bkg_width"])

                # do the convolution
                mode = from_currsys(self.meta["convolve_mode"], self.cmds)

                if image.ndim == 2 and kernel.ndim == 2:
                    new_image = convolve(image - bkg_level, kernel, mode=mode)
                elif image.ndim == 3 and kernel.ndim == 2:
                    kernel = kernel[None, :, :]
                    bkg_level = bkg_level[:, None, None]
                    new_image = convolve(image - bkg_level, kernel, mode=mode)
                elif image.ndim == 3 and kernel.ndim == 3:
                    bkg_level = bkg_level[:, None, None]
                    new_image = np.zeros(image.shape)  # assumes mode="same"
                    for iplane in range(image.shape[0]):
                        new_image[iplane,] = convolve(
                            image[iplane,] - bkg_level[iplane,],
                            kernel[iplane,], mode=mode)

                obj.hdu.data = new_image + bkg_level

                # TODO: careful with which dimensions mean what
                d_x = new_image.shape[-1] - image.shape[-1]
                d_y = new_image.shape[-2] - image.shape[-2]
                for wcsid in ["", "D"]:
                    if "CRPIX1" + wcsid in obj.hdu.header:
                        obj.hdu.header["CRPIX1" + wcsid] += d_x / 2
                        obj.hdu.header["CRPIX2" + wcsid] += d_y / 2

        return obj

    def get_kernel(self, obj):
        self.valid_waverange = None
        if self.kernel is None:
            self.kernel = np.ones((1, 1))
        return self.kernel

    @staticmethod
    def _round_kernel_edges(kernel: np.ndarray) -> np.ndarray:
        x, y = np.array(kernel.shape) // 2
        threshold = min(kernel[x, 0], kernel[x, -1],
                        kernel[0, y], kernel[-1, y])
        kernel[kernel < threshold] = 0.
        # TODO: maybe masked array here?
        return kernel

    def plot(self, obj=None, **kwargs):
        fig, axes = figure_factory()

        kernel = self.get_kernel(obj)
        axes.imshow(kernel, norm="log", origin="lower", **kwargs)

        return fig


def rotational_blur(image, angle):
    """
    Rotate and coadd an image over a given angle to imitate a blur.

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

    n_angles = 1
    rad_to_deg = 57.29578
    edge_pixel_unit_angle = np.arctan2(1, (image.shape[0] // 2)) * rad_to_deg
    while abs(angle) > edge_pixel_unit_angle and n_angles < 25:
        angle /= 2.
        image_rot += rotate(image_rot, angle, reshape=False, order=1)
        # each time kernel is rotated and added, the frame total doubles
        n_angles *= 2

    return image_rot / n_angles


def get_bkg_level(obj, bg_w):
    """
    Determine the background level of image or cube slices.

    Returns a scalar if obj is a 2d image or a vector if obj is a 3D cube (one
    value for each plane).
    The method for background determination is decided by
    self.meta["bkg_width"]:
    If 0, the background is returned as zero (implying no background
    subtraction).
    If -1, the background is estimated as the median of the entire image (or
    cube plane).
    If positive, the background is estimated as the median of a frame of width
    `bkg_width` around the edges.
    """
    if obj.ndim not in (2, 3):
        raise ValueError("Unsupported dimension:", obj.ndim)

    if bg_w == 0:
        if obj.ndim == 3:
            return np.array([0] * obj.shape[0])
        return 0.  # ndim == 2

    mask = np.zeros_like(obj, dtype=bool)
    mask[..., bg_w:-bg_w, bg_w:-bg_w] = True
    bkg = np.ma.masked_array(obj, mask=mask)

    if obj.ndim == 3:
        return np.ma.median(bkg, axis=(2, 1)).data
    return np.ma.median(bkg)  # ndim == 2
