# -*- coding: utf-8 -*-
"""."""

import numpy as np
from scipy.signal import convolve
from astropy import units as u

from ..effects import Effect
from ...base_classes import ImagePlaneBase, FieldOfViewBase, FOVSetupBase
from ...utils import from_currsys, quantify, figure_factory, get_logger
from . import psf_utils as pu


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
    def __init__(self, **kwargs):
        self.kernel = None
        self.valid_waverange = None
        self._waveset = []
        super().__init__(**kwargs)

        params = {
            "flux_accuracy": "!SIM.computing.flux_accuracy",
            "sub_pixel_flag": "!SIM.sub_pixel.flag",
            "z_order": [40, 640],
            "convolve_mode": "same",      # "full", "same"
            "bkg_width": -1,
            "wave_key": "WAVE0",
            "normalise_kernel": True,
            "rotational_blur_angle": 0,
            "report_plot_include": True,
            "report_table_include": False,
        }
        self.meta.update(params)
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta, self.cmds)
        self.convolution_classes = (FieldOfViewBase, ImagePlaneBase)

    def apply_to(self, obj, **kwargs):
        """Apply the PSF."""
        # 1. During setup of the FieldOfViews
        if isinstance(obj, FOVSetupBase) and self._waveset is not None:
            waveset = self._waveset
            if len(waveset) != 0:
                waveset_edges = 0.5 * (waveset[:-1] + waveset[1:])
                obj.split("wave", quantify(waveset_edges, u.um).value)

        # 2. During observe: convolution
        elif isinstance(obj, self.convolution_classes):
            if ((hasattr(obj, "fields") and len(obj.fields) > 0) or
                    (obj.hdu is not None)):
                kernel = self.get_kernel(obj).astype(float)

                # apply rotational blur for field-tracking observations
                rot_blur_angle = self.meta["rotational_blur_angle"]
                if abs(rot_blur_angle) > 0:
                    # makes a copy of kernel
                    kernel = pu.rotational_blur(kernel, rot_blur_angle)

                # normalise psf kernel      KERNEL SHOULD BE normalised within get_kernel()
                # if from_currsys(self.meta["normalise_kernel"], self.cmds):
                #    kernel /= np.sum(kernel)
                #    kernel[kernel < 0.] = 0.

                image = obj.hdu.data.astype(float)

                # subtract background level before convolving, re-add afterwards
                bkg_level = pu.get_bkg_level(image, self.meta["bkg_width"])

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

    def fov_grid(self, which="waveset", **kwargs):
        """See parent docstring."""
        waveset = []
        if which == "waveset":
            if self._waveset is not None:
                _waveset = self._waveset
                waves = 0.5 * (np.array(_waveset)[1:] +
                               np.array(_waveset)[:-1])
                wave_min = kwargs.get("wave_min", np.min(_waveset))
                wave_max = kwargs.get("wave_max", np.max(_waveset))
                mask = (wave_min < waves) * (waves < wave_max)
                waveset = np.unique([wave_min] + list(waves[mask]) +
                                    [wave_max])

        return waveset

    def get_kernel(self, obj):
        self.valid_waverange = None
        if self.kernel is None:
            self.kernel = np.ones((1, 1))
        return self.kernel

    def plot(self, obj=None, **kwargs):
        fig, axes = figure_factory()

        kernel = self.get_kernel(obj)
        axes.imshow(kernel, norm="log", origin="lower", **kwargs)

        return fig
