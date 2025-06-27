# -*- coding: utf-8 -*-
"""Effects for the METIS ifu_cube mode."""

from typing import ClassVar

import numpy as np
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel

from ..effects import Effect
from ...optics.fov import FieldOfView
from ...utils import from_currsys, get_logger

logger = get_logger(__name__)


class LineSpreadFunction(Effect):
    """
    Compute and apply line spread function to IFU cube.

    The effect can be instantiated either with the single parameter

    - lsfwidth : float
         Width of the line-spread function in pixels

    or with the following set of parameters:

    - wavelen : float
         Central wavelength of the LMS setting [um]
    - fit_slope, fit_intercept : float
         Parameters of a linear fit of the mean dispersion (Delta lambda
         per pixel on the LMS detector) as a function of central wavelength.
         For METIS LMS, these values are obtained from `TRACE_LMS.fits` with
         the script `fit_ifu_dispersion.py` in `irdb/METIS/code`.
    - slice_width : float
         On-sky width of a slice of the LMS slicer [arcsec]
    - pixel_scale : float
         On-sky angle covered by a pixel of the LMS detector array [arcsec]
    - spec_binwidth : float
         Spectral bin width of the 3D detector of the lms_cube mode.

    These values are set in `METIS_LMS_SMPL.yaml`.

    .. versionadded:: 0.10.0

    """

    z_order: ClassVar[tuple[int, ...]] = (660,)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta, self.cmds)

        if "lsfwidth" not in self.meta:
            self.lsfwidth = self.get_lsf_width()
        else:
            self.lsfwidth = self.meta["lsfwidth"]
        self.kernel = self.get_kernel()

    def apply_to(self, obj, **kwargs):
        """Apply the LSF."""
        if not isinstance(obj, FieldOfView):
            return obj

        if obj.hdu is None or obj.hdu.header["NAXIS"] != 3:
            logger.error("Cannot apply LSF convolution")
            return obj

        obj.hdu.data = convolve(obj.hdu.data, self.kernel)
        obj.cube = obj.hdu
        return obj

    def get_lsf_width(self):
        """Determine width of the LSF kernel at central wavelength."""
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
        """Build LSF kernel: box kernel smoothed with narrow Gauss."""
        box = Box1DKernel(width=self.lsfwidth)
        gauss = Gaussian1DKernel(1)
        if box.shape > gauss.shape:
            kernel = convolve(box.array, gauss)[:, np.newaxis, np.newaxis]
        else:
            kernel = convolve(gauss.array, box)[:, np.newaxis, np.newaxis]
        return kernel
