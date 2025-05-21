# -*- coding: utf-8 -*-
"""Effects for the METIS ifu_cube mode"""

from typing import ClassVar

import numpy as np
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel

from ..effects import Effect
from ...optics.fov import FieldOfView
from ...utils import from_currsys, get_logger

logger = get_logger(__name__)


class LineSpreadFunction(Effect):
    """
    Compute and apply line spread function to IFU cube
    """
    z_order: ClassVar[tuple[int, ...]] = (660,)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {

        }
        self.meta.update(params)
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta, self.cmds)

        if "lsfwidth" not in self.meta:
            self.lsfwidth = self.get_lsf_width()
        else:
            self.lsfwidth = self.meta["lsfwidth"]
        self.kernel = self.get_kernel()

    def apply_to(self, fov, **kwargs):
        """Apply the LSF"""
        if not isinstance(fov, FieldOfView):
            return fov

        if fov.hdu is None or fov.hdu.header["NAXIS"] != 3:
            logger.error("Cannot apply LSF convolution")
            return fov


        print("Happily convolving with LSF")
        fov.hdu.data = convolve(fov.hdu.data, self.kernel)
        fov.cube = fov.hdu
        return fov

    def get_lsf_width(self):
        """Determine width of the LSF kernel at central wavelength"""

        slope = self.meta['fit_slope']
        intercept = self.meta['fit_intercept']
        lamc = self.meta["wavelen"]
        dlam_per_pix = slope * lamc + intercept

        slice_width = self.meta['slice_width']
        pixel_scale = self.meta['pixel_scale']

        dlam_per_slice = dlam_per_pix * slice_width / pixel_scale

        spec_binwidth = self.meta["spec_binwidth"]

        return dlam_per_slice / spec_binwidth


    def get_kernel(self):
        """Build LSF kernel: box kernel smoothed with narrow Gauss"""

        box = Box1DKernel(width=self.lsfwidth)
        gauss = Gaussian1DKernel(1)
        if box.shape > gauss.shape:
            kernel = convolve(box.array, gauss)[np.newaxis, np.newaxis, :]
        else:
            kernel = convolve(gauss.array, box)[np.newaxis, np.newaxis, :]
        return kernel
