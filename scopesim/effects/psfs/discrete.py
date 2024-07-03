# -*- coding: utf-8 -*-
"""."""

import numpy as np
from scipy.signal import convolve
from scipy.interpolate import RectBivariateSpline

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS


from ...base_classes import FieldOfViewBase
from ...utils import from_currsys, check_keys
from . import PSF, PoorMansFOV
from . import psf_utils as pu


class DiscretePSF(PSF):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [43]
        self.convolution_classes = FieldOfViewBase
        # self.convolution_classes = ImagePlaneBase


class FieldConstantPSF(DiscretePSF):
    """A PSF that is constant across the field.

    For spectroscopy, a wavelength-dependent PSF cube is built, where for each
    wavelength the reference PSF is scaled proportional to wavelength.
    """

    required_keys = {"filename"}

    def __init__(self, **kwargs):
        # sub_pixel_flag and flux_accuracy are taken care of in PSF base class
        super().__init__(**kwargs)

        check_keys(self.meta, self.required_keys, action="error")

        self.meta["z_order"] = [262, 662]
        self._waveset, self.kernel_indexes = pu.get_psf_wave_exts(
            self._file, self.meta["wave_key"])
        self.current_layer_id = None
        self.current_ext = None
        self.current_data = None
        self.kernel = None

    def get_kernel(self, fov):
        """Find nearest wavelength and build PSF kernel from file"""
        idx = pu.nearest_index(fov.wavelength, self._waveset)
        ext = self.kernel_indexes[idx]
        if ext != self.current_layer_id:
            if fov.hdu.header["NAXIS"] == 3:
                self.current_layer_id = ext
                self.make_psf_cube(fov)
            else:
                self.kernel = self._file[ext].data
                self.current_layer_id = ext
                hdr = self._file[ext].header

                self.kernel /= np.sum(self.kernel)

                # compare kernel and fov pixel scales, rescale if needed
                if "CUNIT1" in hdr:
                    unit_factor = u.Unit(hdr["CUNIT1"].lower()).to(u.deg)
                else:
                    unit_factor = 1

                kernel_pixel_scale = hdr["CDELT1"] * unit_factor
                fov_pixel_scale = fov.header["CDELT1"]

                # rescaling kept inside loop to avoid rescaling for every fov
                pix_ratio = kernel_pixel_scale / fov_pixel_scale
                if abs(pix_ratio - 1) > self.meta["flux_accuracy"]:
                    spline_order = from_currsys(
                        "!SIM.computing.spline_order", cmds=self.cmds)
                    self.kernel = pu.rescale_kernel(
                        self.kernel, pix_ratio, spline_order)

                if ((fov.header["NAXIS1"] < hdr["NAXIS1"]) or
                        (fov.header["NAXIS2"] < hdr["NAXIS2"])):
                    self.kernel = pu.cutout_kernel(self.kernel, fov.header,
                                                   kernel_header=hdr)

        return self.kernel

    def make_psf_cube(self, fov):
        """Create a wavelength-dependent psf cube"""

        # Some data from the fov
        nxfov, nyfov = fov.hdu.header["NAXIS1"], fov.hdu.header["NAXIS2"]
        fov_pixel_scale = fov.hdu.header["CDELT1"]
        fov_pixel_unit = fov.hdu.header["CUNIT1"].lower()

        lam = fov.hdu.header["CDELT3"] * (1 + np.arange(fov.hdu.header["NAXIS3"])
                                          - fov.hdu.header["CRPIX3"]) \
                                          + fov.hdu.header["CRVAL3"]

        # adapt the size of the output cube to the FOV's spatial shape
        nxpsf = min(512, 2 * nxfov + 1)
        nypsf = min(512, 2 * nyfov + 1)

        # Some data from the psf file
        ext = self.current_layer_id
        hdr = self._file[ext].header
        refwave = hdr[self.meta["wave_key"]]

        if "CUNIT1" in hdr:
            unit_factor = u.Unit(hdr["CUNIT1"].lower()).to(
                u.Unit(fov_pixel_unit))
        else:
            unit_factor = 1
        ref_pixel_scale = hdr["CDELT1"] * unit_factor

        psfwcs = WCS(hdr)
        psf = self._file[ext].data
        psf = psf/psf.sum()         # normalisation of the input psf
        nxin, nyin = psf.shape

        # We need linear interpolation to preserve positivity. Might think of
        # more elaborate positivity-preserving schemes.
        # Note: According to some basic profiling, this line is one of the
        #       single largest hits on performance.
        ipsf = RectBivariateSpline(np.arange(nyin), np.arange(nxin), psf,
                                   kx=1, ky=1)

        xcube, ycube = np.meshgrid(np.arange(nxpsf), np.arange(nypsf))
        cubewcs = WCS(naxis=2)
        cubewcs.wcs.ctype = ["LINEAR", "LINEAR"]
        cubewcs.wcs.crval = [0., 0.]
        cubewcs.wcs.crpix = [(nxpsf + 1) / 2, (nypsf + 1) / 2]
        cubewcs.wcs.cdelt = [fov_pixel_scale, fov_pixel_scale]
        cubewcs.wcs.cunit = [fov_pixel_unit, fov_pixel_unit]

        xworld, yworld = cubewcs.all_pix2world(xcube, ycube, 1)
        outcube = np.zeros((lam.shape[0], nypsf, nxpsf), dtype=np.float32)
        for i, wave in enumerate(lam):
            psf_wave_pixscale = ref_pixel_scale * wave / refwave
            psfwcs.wcs.cdelt = [psf_wave_pixscale,
                                psf_wave_pixscale]
            xpsf, ypsf = psfwcs.all_world2pix(xworld, yworld, 0)
            outcube[i,] = (ipsf(ypsf, xpsf, grid=False)
                           * fov_pixel_scale**2 / psf_wave_pixscale**2)

        self.kernel = outcube.reshape((lam.shape[0], nypsf, nxpsf))
        # fits.writeto("test_psfcube.fits", data=self.kernel, overwrite=True)

    def plot(self):
        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        spec_dict = from_currsys("!SIM.spectral", self.cmds)
        return super().plot(PoorMansFOV(pixel_scale, spec_dict))


class FieldVaryingPSF(DiscretePSF):
    """
    TBA.

    Parameters
    ----------
    sub_pixel_flag : bool, optional
    flux_accuracy : float, optional
        Default 1e-3. Level of flux conservation during rescaling of kernel

    """

    required_keys = {"filename"}

    def __init__(self, **kwargs):
        # sub_pixel_flag and flux_accuracy are taken care of in PSF base class
        super().__init__(**kwargs)

        check_keys(self.meta, self.required_keys, action="error")

        self.meta["z_order"] = [261, 661]

        ws, ki = pu.get_psf_wave_exts(self._file, self.meta["wave_key"])
        self._waveset, self.kernel_indexes = ws, ki
        self.current_ext = None
        self.current_data = None
        self._strehl_imagehdu = None

    def apply_to(self, fov, **kwargs):
        """See parent docstring."""
        # TODO: add in field rotation
        # TODO: add in 3D cubes
        # accept "full", "dit", "none"

        # check if there are any fov.fields to apply a psf to
        if isinstance(fov, FieldOfViewBase):
            if len(fov.fields) > 0:
                if fov.image is None:
                    fov.image = fov.make_image_hdu().data

                old_shape = fov.image.shape

                # Get kernels that cover this fov, and their respective masks.
                # Kernels and masks returned by .get_kernel as list of tuples.
                canvas = None
                kernels_masks = self.get_kernel(fov)
                for kernel, mask in kernels_masks:

                    # renormalise the kernel if needs be
                    kernel[kernel < 0.] = 0.
                    sum_kernel = np.sum(kernel)
                    if abs(sum_kernel - 1) > self.meta["flux_accuracy"]:
                        kernel /= sum_kernel

                    # image convolution
                    image = fov.image.astype(float)
                    kernel = kernel.astype(float)
                    new_image = convolve(image, kernel, mode="same")
                    if canvas is None:
                        canvas = np.zeros(new_image.shape)

                    # mask convolution + combine with convolved image
                    if mask is not None:
                        new_mask = convolve(mask, kernel, mode="same")
                        canvas += new_image * new_mask
                    else:
                        canvas = new_image

                # reset WCS header info
                new_shape = canvas.shape
                fov.image = canvas

                # TODO: careful with which dimensions mean what
                if "CRPIX1" in fov.header:
                    fov.header["CRPIX1"] += (new_shape[0] - old_shape[0]) / 2
                    fov.header["CRPIX2"] += (new_shape[1] - old_shape[1]) / 2

                if "CRPIX1D" in fov.header:
                    fov.header["CRPIX1D"] += (new_shape[0] - old_shape[0]) / 2
                    fov.header["CRPIX2D"] += (new_shape[1] - old_shape[1]) / 2

        return fov

    def get_kernel(self, fov):
        # 0. get file extension
        # 1. pull out strehl map for fov header
        # 2. get number of unique psfs
        # 3. pull out those psfs
        # 4. if more than one, make masks for the fov on the fov pixel scale
        # 5. make list of tuples with kernel and mask

        # find which file extension to use - keep pointer in self.current_data
        fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        jj = pu.nearest_index(fov_wave, self._waveset)
        ext = self.kernel_indexes[jj]
        if ext != self.current_ext:
            self.current_ext = ext
            self.current_data = self._file[ext].data

        # compare the fov and psf pixel scales
        kernel_pixel_scale = self._file[ext].header["CDELT1"]
        fov_pixel_scale = fov.header["CDELT1"]

        # get the spatial map of the kernel cube layers
        strl_hdu = self.strehl_imagehdu
        strl_cutout = pu.get_strehl_cutout(fov.header, strl_hdu)

        # get the kernels and mask that fit inside the fov boundaries
        layer_ids = np.round(np.unique(strl_cutout.data)).astype(int)
        if len(layer_ids) > 1:
            kernels = [self.current_data[ii] for ii in layer_ids]
            masks = [strl_cutout.data == ii for ii in layer_ids]
            self.kernel = [[krnl, msk] for krnl, msk in zip(kernels, masks)]
        else:
            self.kernel = [[self.current_data[layer_ids[0]], None]]

        # TODO: re-scale kernel and masks to pixel_scale of FOV
        # TODO: can this be put somewhere else to save on iterations?
        # TODO: should the mask also be rescaled?
        # rescale the pixel scale of the kernel to match the fov images
        pix_ratio = fov_pixel_scale / kernel_pixel_scale
        if abs(pix_ratio - 1) > self.meta["flux_accuracy"]:
            spline_order = from_currsys(
                "!SIM.computing.spline_order", cmds=self.cmds)
            for ii, kern in enumerate(self.kernel):
                self.kernel[ii][0] = pu.rescale_kernel(
                    kern[0], pix_ratio, spline_order)

        for i, kern in enumerate(self.kernel):
            self.kernel[i][0] /= np.sum(kern[0])

        return self.kernel

    @property
    def strehl_imagehdu(self):
        """The HDU containing the positional info for kernel layers."""
        if self._strehl_imagehdu is None:
            ecat = self._file[0].header["ECAT"]
            if isinstance(self._file[ecat], fits.ImageHDU):
                self._strehl_imagehdu = self._file[ecat]

            # TODO: impliment this case
            elif isinstance(self._file[ecat], fits.BinTableHDU):
                cat = self._file[ecat]
                self._strehl_imagehdu = pu.make_strehl_map_from_table(cat)

        return self._strehl_imagehdu

    def plot(self):
        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        spec_dict = from_currsys("!SIM.spectral", self.cmds)
        return super().plot(PoorMansFOV(pixel_scale, spec_dict))
