# -*- coding: utf-8 -*-
"""Contains the base class for all PSF effects."""

from typing import ClassVar

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import convolve
from scipy.ndimage import rotate
from astropy import units as u

from ..effects import Effect
from ...optics import ImagePlane
from ...optics.fov import FieldOfView
from ...optics.fov_volume_list import FovVolumeList
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
            "rotational_blur_angle": 0*u.deg,
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        self.meta = from_currsys(self.meta, self.cmds)
        self.convolution_classes = (FieldOfView, ImagePlane)

    def __call__(self, data: ArrayLike, kernel: ArrayLike) -> NDArray:
        # subtract background level before convolving, re-add afterwards
        bkg_level = get_bkg_level(data, self.meta["bkg_width"])

        if data.ndim == 3:
            bkg_level = bkg_level[:, None, None]

        image = np.asarray(data - bkg_level)

        # Shortcut for flat images, which wouldn't change with convolution
        if image.sum() == 0.:
            logger.debug(
                "%s: uniform field, convolution skipped", self.display_name)
            return data  # return unchanged

        mode = from_currsys(self.meta["convolve_mode"], self.cmds)

        match image.ndim, kernel.ndim:
            case 2, 2:
                return convolve(image, kernel, mode=mode) + bkg_level

            case 3, 2:
                kernel = kernel[None, ...]
                return convolve(image, kernel, mode=mode) + bkg_level

            case 3, 3:
                if mode != "same":
                    logger.warning(
                        "cube-cube convolution assumes mode='same', but mode "
                        "is %s, result may be inconsistent.", mode)

                convolved = np.zeros(data.shape)  # assumes mode="same"
                for iplane in range(data.shape[0]):
                    convolved[iplane,] = convolve(
                        image[iplane,],
                        kernel[iplane,],
                        mode=mode,
                    )
                return convolved + bkg_level

            case _:
                raise ValueError("Image and Kernel shape mismatch.")

    def _apply_to_fvl(self, fvl: FovVolumeList) -> None:
        logger.debug("Apply %s, FoV setup", self.display_name)

        if len(self._waveset) == 0:
            return

        waveset_edges = 0.5 * (self._waveset[:-1] + self._waveset[1:])
        fvl.split("wave", quantify(waveset_edges, u.um).value)

    def _apply_convolution(self, obj) -> None:
        logger.debug("Apply %s, convolution", self.display_name)

        if ((not hasattr(obj, "fields") or len(obj.fields) == 0)
            and (obj.hdu is None)):
            return

        kernel = self.get_kernel(obj).astype(float)

        # This doesn't work because of a "Delta PSF" in some mocks...
        # if kernel.size == 1:  # only 1 pixel
        #     raise ValueError("Cannot convolve single pixel PSF.")

        # apply rotational blur for field-tracking observations
        rot_blur_angle = self.meta["rotational_blur_angle"]
        if abs(rot_blur_angle << u.deg) > 0*u.deg:
            # makes a copy of kernel
            kernel = rotational_blur(kernel, rot_blur_angle)

        # Round the edges of kernels to avoid square stars
        if self.meta.get("rounded_edges", False) and kernel.ndim == 2:
            kernel = self._round_kernel_edges(kernel)

        # normalise psf kernel      KERNEL SHOULD BE normalised within get_kernel()
        # if from_currsys(self.meta["normalise_kernel"], self.cmds):
        #    kernel /= np.sum(kernel)
        #    kernel[kernel < 0.] = 0.

        orig_shape = obj.hdu.data.shape

        logger.debug("PSF convolution start")
        obj.hdu.data = self(obj.hdu.data.astype(float), kernel)
        logger.debug("PSF convolution done")

        # TODO: careful with which dimensions mean what
        d_x = obj.hdu.data.shape[-1] - orig_shape[-1]
        d_y = obj.hdu.data.shape[-2] - orig_shape[-2]
        for wcsid in ["", "D"]:
            if "CRPIX1" + wcsid in obj.hdu.header:
                obj.hdu.header["CRPIX1" + wcsid] += d_x / 2
                obj.hdu.header["CRPIX2" + wcsid] += d_y / 2

    def apply_to(self, obj, **kwargs):
        """Apply the PSF."""
        # 1. During setup of the FieldOfViews
        if isinstance(obj, FovVolumeList) and self._waveset is not None:
            self._apply_to_fvl(obj)

        # 2. During observe: convolution
        elif isinstance(obj, self.convolution_classes):
            self._apply_convolution(obj)

        return obj

    def get_kernel(self, obj):
        self.valid_waverange = None
        if self.kernel is None:
            self.kernel = np.ones((1, 1))
        return self.kernel

    @staticmethod
    def _round_kernel_edges(kernel: ArrayLike) -> NDArray:
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


@u.quantity_input
def rotational_blur(image, angle: u.Quantity[u.deg]):
    """
    Rotate and coadd an image over a given angle to imitate a blur.

    Parameters
    ----------
    image : array-like
        Image to blur.
    angle : u.Quantity["angle"]
        Angle over which the image should be rotationally blurred.

    .. versionchanged:: 0.11.3

       Require `angle` to be a Quantity.

    Returns
    -------
    image_rot : NDArray
        Blurred image

    """
    image_rot = np.copy(image)

    edge_pixel_unit_angle = np.arctan2(1, (image.shape[0] // 2)) * u.rad
    n_steps = np.ceil(np.log2(abs(angle) / edge_pixel_unit_angle))
    n_steps = int(min(n_steps, 8))  # avoid overrun

    current_angle = angle.copy()
    for _ in range(n_steps):
        current_angle /= 2.
        image_rot += rotate(image_rot, current_angle, reshape=False, order=3)
        # each time kernel is rotated and added, the frame total doubles

    return image_rot / image_rot.sum() * image.sum()


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
