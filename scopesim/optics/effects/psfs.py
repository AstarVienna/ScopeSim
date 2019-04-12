from copy import deepcopy

import numpy as np
from scipy.signal import convolve
from scipy.ndimage import zoom
from scipy.interpolate import griddata

from astropy.io import fits
from astropy import units as u

from ... import utils
from ... import rc
from ..fov import FieldOfView
from .. import image_plane_utils as imp_utils
from .effects import Effect


class PSF(Effect):
    def __init__(self, **kwargs):
        self.kernel = None
        self.valid_waverange = None
        super(PSF, self).__init__(**kwargs)
        self.meta["SIM_FLUX_ACCURACY"] = rc.__rc__["SIM_FLUX_ACCURACY"]
        self.meta["SIM_SUB_PIXEL_ACCURACY"] = rc.__rc__["SIM_SUB_PIXEL_ACCURACY"]

        self.meta.update(kwargs)

    def apply_to(self, obj):
        if len(obj.fields) > 0:
            kernel = self.get_kernel(obj)

            sub_pixel = False
            if "SIM_SUB_PIXEL_ACCURACY" in self.meta:
                sub_pixel = self.meta["SIM_SUB_PIXEL_ACCURACY"]

            if obj.hdu.data is None:
                obj.view(sub_pixel)

            old_shape = obj.hdu.data.shape
            new_image = convolve(obj.hdu.data, kernel, mode="full")
            new_shape = new_image.shape

            obj.hdu.data = new_image

            # ..todo: careful with which dimensions mean what
            if "CRPIX1" in obj.hdu.header:
                obj.hdu.header["CRPIX1"] += (new_shape[0] - old_shape[0]) / 2
                obj.hdu.header["CRPIX2"] += (new_shape[1] - old_shape[1]) / 2

            if "CRPIX1D" in obj.hdu.header:
                obj.hdu.header["CRPIX1D"] += (new_shape[0] - old_shape[0]) / 2
                obj.hdu.header["CRPIX2D"] += (new_shape[1] - old_shape[1]) / 2

        return obj

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"wavelengths": waverange}

    def get_kernel(self, fov):
        self.valid_waverange = None
        self.kernel = np.ones((1, 1))
        return self.kernel


################################################################################
# Analytical PSFs - Vibration, Seeing, NCPAs


class AnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(AnalyticalPSF, self).__init__(**kwargs)


class Vibration(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(Vibration, self).__init__(**kwargs)
        self.meta["z_order"] = [400]


class NonCommonPathAberration(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(NonCommonPathAberration, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]


class Seeing(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(Seeing, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]


################################################################################
# Semi-analytical PSFs - Poppy PSFs


class SemiAnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(SemiAnalyticalPSF, self).__init__(**kwargs)


class PoppyFieldVaryingPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(PoppyFieldVaryingPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]


class PoppyFieldConstantPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(PoppyFieldConstantPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]


################################################################################
# Discreet PSFs - MAORY and co PSFs


class DiscretePSF(PSF):
    def __init__(self, **kwargs):
        super(DiscretePSF, self).__init__(**kwargs)


class FieldConstantPSF(DiscretePSF):
    def __init__(self, **kwargs):
        super(FieldConstantPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]
        self.waveset, self.kernel_indexes = get_psf_wave_exts(self)
        self.current_layer_id = None

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"wavelengths": self.waveset}

    def get_kernel(self, fov):
        # fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        # ii = nearest_index(fov_wave, self.waveset)
        # ext = self.kernel_indexes[ii]
        # if ext != self.current_layer_id:
        #     self.kernel = self._file[ext]
        #     self.current_layer_id = ext
        self.kernel = self._file[2].data[0]


class FieldVaryingPSF(DiscretePSF):
    def __init__(self, **kwargs):
        super(FieldVaryingPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]
        self.waveset, self.kernel_indexes = get_psf_wave_exts(self._file)
        self.current_ext = None
        self.current_data = None
        self._strehl_imagehdu = None

    def apply_to(self, fov):
        if len(fov.fields) > 0:
            if fov.hdu.data is None:
                fov.view(self.meta["SIM_SUB_PIXEL_ACCURACY"])

            old_shape = fov.hdu.data.shape

            canvas = None
            kernels_masks = self.get_kernel(fov)
            for kernel, mask in kernels_masks:

                sum_kernel = np.sum(kernel)
                if abs(sum_kernel - 1) > self.meta["SIM_FLUX_ACCURACY"]:
                    kernel /= sum_kernel

                new_image = convolve(fov.hdu.data, kernel, mode="same")
                if canvas is None:
                    canvas = np.zeros(new_image.shape)

                if mask is not None:
                    new_mask =  convolve(mask, kernel, mode="same")
                    canvas += new_image * new_mask
                else:
                    canvas = new_image

            new_shape = canvas.shape
            fov.hdu.data = canvas

            # ..todo: careful with which dimensions mean what
            if "CRPIX1" in fov.hdu.header:
                fov.hdu.header["CRPIX1"] += (new_shape[0] - old_shape[0]) / 2
                fov.hdu.header["CRPIX2"] += (new_shape[1] - old_shape[1]) / 2

            if "CRPIX1D" in fov.hdu.header:
                fov.hdu.header["CRPIX1D"] += (new_shape[0] - old_shape[0]) / 2
                fov.hdu.header["CRPIX2D"] += (new_shape[1] - old_shape[1]) / 2

        return fov

    def get_kernel(self, fov):
        # 0. get file extension
        # 1. pull out strehl map for fov header
        # 2. get number of unique psfs
        # 3. pull out those psfs
        # 4. if more than one, make masks for the fov on the fov pixel scale
        # 5. make list of tuples with kernel and mask

        fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        jj = nearest_index(fov_wave, self.waveset)
        ii = self.kernel_indexes[jj]
        if ii != self.current_ext:
            self.current_ext = ii
            self.current_data = self._file[ii].data
        kernel_pixel_scale = self._file[ii].header["CDELT1"]
        fov_pixel_scale = fov.hdu.header["CDELT1"]

        strl_hdu = self.strehl_imagehdu
        strl_cutout = get_strehl_cutout(fov.hdu.header, strl_hdu)

        layer_ids = np.round(np.unique(strl_cutout.data)).astype(int)
        if len(layer_ids) > 1:
            kernels = [self.current_data[ii] for ii in layer_ids]
            masks = [strl_cutout.data.T == ii for ii in layer_ids]          # there's a .T in here that I don't like
            self.kernel = [[krnl, msk] for krnl, msk in zip(kernels, masks)]
        else:
            self.kernel = [[self.current_data[layer_ids[0]], None]]

        # .. todo: re-scale kernel and masks to pixel_scale of FOV
        # .. todo: can this be put somewhere else to save on iterations?
        pix_ratio = fov_pixel_scale / kernel_pixel_scale
        if abs(pix_ratio - 1) > self.meta["SIM_FLUX_ACCURACY"]:
            for ii in range(len(self.kernel)):
                self.kernel[ii][0] = resize_array(self.kernel[ii][0], pix_ratio)

        return self.kernel

    @property
    def strehl_imagehdu(self):
        if self._strehl_imagehdu is None:
            ecat = self._file[0].header["ECAT"]
            if isinstance(self._file[ecat], fits.ImageHDU):
                self._strehl_imagehdu = self._file[ecat]

            # ..todo: impliment this case
            elif isinstance(self._file[ecat], fits.BinTableHDU):
                cat = self._file[ecat]
                self._strehl_imagehdu = make_strehl_map_from_table(cat)

        return self._strehl_imagehdu


################################################################################
# Helper functions

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




def resize_array(image, scale_factor, order=1):
    sum_image = np.sum(image)
    image = zoom(image, scale_factor, order=order)
    image = np.nan_to_num(image, copy=False)        # numpy version >=1.13
    sum_new_image = np.sum(image)
    image *= sum_image / sum_new_image

    return image


def get_strehl_cutout(fov_header, strehl_imagehdu):

    image = np.zeros((fov_header["NAXIS1"], fov_header["NAXIS2"]))
    canvas_hdu = fits.ImageHDU(header=fov_header, data=image)
    canvas_hdu = imp_utils.add_imagehdu_to_imagehdu(strehl_imagehdu,
                                                    canvas_hdu, order=0,
                                                    conserve_flux=False)
    canvas_hdu.data = canvas_hdu.data.astype(int)

    return canvas_hdu


def nearest_index(x, x_array):
    # return int(round(np.interp(x, x_array, np.arange(len(x_array)))))
    return np.argmin(abs(x_array - x))


def get_psf_wave_exts(hdu_list):
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

    wave_ext = [ii for ii in range(len(hdu_list))
                if "WAVE0" in hdu_list[ii].header]
    wave_set = [hdu.header["WAVE0"] for hdu in hdu_list
                if "WAVE0" in hdu.header]
    wave_set = utils.quantify(wave_set, u.um)

    return wave_set, wave_ext


