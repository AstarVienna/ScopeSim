import numpy as np
from scipy.signal import convolve
from scipy.ndimage import zoom
from scipy.interpolate import griddata

from astropy import units as u
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel

from .effects import Effect
from ..optics import image_plane_utils as imp_utils
from ..base_classes import ImagePlaneBase, FieldOfViewBase
from .. import utils
from .. import rc


class PSF(Effect):
    def __init__(self, **kwargs):
        self.kernel = None
        self.valid_waverange = None
        self._waveset = None
        flux_accuracy = rc.__config__["!SIM.computing.flux_accuracy"]
        super(PSF, self).__init__(**kwargs)

        self.meta["flux_accuracy"] = float(flux_accuracy)
        self.meta["sub_pixel_flag"] = rc.__config__["!SIM.sub_pixel.flag"]
        self.meta.update(kwargs)
        self.apply_to_classes = (FieldOfViewBase, ImagePlaneBase)

    def apply_to(self, obj):
        if isinstance(obj, self.apply_to_classes):
            if (hasattr(obj, "fields") and len(obj.fields) > 0) or \
                    obj.hdu.data is not None:

                if obj.hdu.data is None:
                    obj.view(self.meta["sub_pixel_flag"])
                old_shape = obj.hdu.data.shape

                kernel = self.get_kernel(obj)
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

    def fov_grid(self, which="waveset", **kwargs):
        waves = []
        if which == "waveset":
            waves = self._waveset

        return waves

    def get_kernel(self, obj):
        self.valid_waverange = None
        self.kernel = np.ones((1, 1))
        return self.kernel


################################################################################
# Analytical PSFs - Vibration, Seeing, NCPAs


class AnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(AnalyticalPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [40, 340]


class Vibration(AnalyticalPSF):
    """
    Creates a wavelength independent kernel image
    """
    def __init__(self, **kwargs):
        super(Vibration, self).__init__(**kwargs)
        self.meta["z_order"] = [44, 444]
        self.meta["width_n_fwhms"] = 4
        self.apply_to_classes = ImagePlaneBase

        self.required_keys = ["fwhm", "pixel_scale"]
        utils.check_keys(self.meta, self.required_keys, action="error")
        self.kernel = None

    def get_kernel(self, implane):
        if self.kernel is None:
            utils.from_currsys(self.meta)
            fwhm_pix = self.meta["fwhm"] / self.meta["pixel_scale"]
            sigma = fwhm_pix / 2.35
            width = int(fwhm_pix * self.meta["width_n_fwhms"])
            self.kernel = Gaussian2DKernel(sigma, x_size=width, y_size=width,
                                           mode="center").array

        return self.kernel


class NonCommonPathAberration(AnalyticalPSF):
    """
    Needed: pixel_scale
    Accepted: kernel_width, strehl_drift
    """
    def __init__(self, **kwargs):
        super(NonCommonPathAberration, self).__init__(**kwargs)
        self.meta["z_order"] = [41, 341]
        self.meta["kernel_width"] = None
        self.meta["strehl_drift"] = 0.02
        self.meta["wave_min"] = "!SIM.spectral.wave_min"
        self.meta["wave_max"] = "!SIM.spectral.wave_max"
        self.apply_to_classes = FieldOfViewBase

        self._total_wfe = None

        self.valid_waverange = [0.1 * u.um, 0.2 * u.um]

        self.required_keys = ["pixel_scale"]
        utils.check_keys(self.meta, self.required_keys, action="error")

    def fov_grid(self, which="waveset", **kwargs):

        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = utils.from_currsys(self.meta)

            min_sr = wfe2strehl(self.total_wfe, self.meta["wave_min"])
            max_sr = wfe2strehl(self.total_wfe, self.meta["wave_max"])

            srs = np.arange(min_sr, max_sr, self.meta["strehl_drift"])
            waves = 6.2831853 * self.total_wfe * (-np.log(srs))**-0.5
            waves = utils.quantify(waves, u.um).to(u.um).value
            waves = (list(waves) + [self.meta["wave_max"]]) * u.um
        else:
            waves = [] * u.um

        return waves

    def get_kernel(self, obj):
        waves = obj.meta["wave_min"], obj.meta["wave_max"]

        old_waves = self.valid_waverange
        wave_mid_old = 0.5 * (old_waves[0] + old_waves[1])
        wave_mid_new = 0.5 * (waves[0] + waves[1])
        strehl_old = wfe2strehl(wfe=self.total_wfe, wave=wave_mid_old)
        strehl_new = wfe2strehl(wfe=self.total_wfe, wave=wave_mid_new)

        if np.abs(1 - strehl_old / strehl_new) > self.meta["strehl_drift"]:
            self.valid_waverange = waves
            self.kernel = wfe2gauss(wfe=self.total_wfe, wave=wave_mid_new,
                                    width=self.meta["kernel_width"])
        return self.kernel

    @property
    def total_wfe(self):
        if self._total_wfe is None:
            if self.table is not None:
                self._total_wfe = get_total_wfe_from_table(self.table)
            else:
                self._total_wfe = 0
        return self._total_wfe


class Seeing(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(Seeing, self).__init__(**kwargs)
        self.meta["z_order"] = [43, 343]


class GaussianDiffractionPSF(AnalyticalPSF):
    def __init__(self, diameter, **kwargs):
        super(GaussianDiffractionPSF, self).__init__(**kwargs)
        self.meta["diameter"] = diameter
        self.meta["z_order"] = [42, 342]

    def fov_grid(self, which="waveset", **kwargs):
        wavelengths = []
        if which == "waveset" and \
                "waverange" in kwargs and \
                "pixel_scale" in kwargs:
            waverange = utils.quantify(kwargs["waverange"], u.um)
            diameter = utils.quantify(self.meta["diameter"], u.m).to(u.um)
            fwhm = 1.22 * (waverange / diameter).value  # in rad

            pixel_scale = utils.quantify(kwargs["pixel_scale"], u.deg)
            pixel_scale = pixel_scale.to(u.rad).value
            fwhm_range = np.arange(fwhm[0], fwhm[1], pixel_scale)
            wavelengths = list(fwhm_range / 1.22 * diameter.to(u.m))

        # .. todo: check that this is actually correct
        return wavelengths

    def update(self, **kwargs):
        if "diameter" in kwargs:
            self.meta["diameter"] = kwargs["diameter"]

    def get_kernel(self, fov):

        pixel_scale = fov.header["CDELT1"] * u.deg.to(u.arcsec)
        pixel_scale = utils.quantify(pixel_scale, u.arcsec)

        wave = 0.5 * (fov.meta["wave_max"] + fov.meta["wave_min"])

        wave = utils.quantify(wave, u.um)
        diameter = utils.quantify(self.meta["diameter"], u.m).to(u.um)
        fwhm = 1.22 * (wave / diameter) * u.rad.to(u.arcsec) / pixel_scale

        sigma = fwhm.value / 2.35
        kernel = Gaussian2DKernel(sigma, mode="center").array

        return kernel


################################################################################
# Semi-analytical PSFs - Poppy PSFs


class SemiAnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(SemiAnalyticalPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [50, 350]


class PoppyFieldVaryingPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(PoppyFieldVaryingPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [51, 351]


class PoppyFieldConstantPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(PoppyFieldConstantPSF, self).__init__(**kwargs)
        self.meta["z_order"] = [52, 352]


################################################################################
# Discreet PSFs - MAORY and co PSFs


class DiscretePSF(PSF):
    def __init__(self, **kwargs):
        super(DiscretePSF, self).__init__(**kwargs)
        self.meta["z_order"] = [60, 360]


class FieldConstantPSF(DiscretePSF):
    def __init__(self, **kwargs):
        # sub_pixel_flag and flux_accuracy are taken care of in PSF base class
        super(FieldConstantPSF, self).__init__(**kwargs)

        self.required_keys = ["filename"]
        utils.check_keys(self.meta, self.required_keys, action="error")

        self.meta["z_order"] = [62, 362]
        self._waveset, self.kernel_indexes = get_psf_wave_exts(self._file)
        self.current_layer_id = None
        self.current_ext = None
        self.current_data = None
        self.kernel = None

    # def apply_to(self, fov):
    #   Taken care of by PSF base class

    # def fov_grid(self, which="waveset", **kwargs):
    #   Taken care of by PSF base class

    def get_kernel(self, fov):
        # find nearest wavelength and pull kernel from file
        fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        ii = nearest_index(fov_wave, self._waveset)
        ext = self.kernel_indexes[ii]
        if ext != self.current_layer_id:
            self.kernel = self._file[ext].data
            self.current_layer_id = ext

            # compare kernel and fov pixel scales, rescale if needed
            kernel_pixel_scale = self._file[ext].header["CDELT1"]
            fov_pixel_scale = fov.hdu.header["CDELT1"]

            # rescaling kept inside loop to avoid rescaling for every fov
            pix_ratio = kernel_pixel_scale / fov_pixel_scale
            if abs(pix_ratio - 1) > self.meta["flux_accuracy"]:
                self.kernel = resize_array(self.kernel, pix_ratio)

        return self.kernel


class FieldVaryingPSF(DiscretePSF):
    """
    kwargs
    ------
    sub_pixel_flag
    flux_accuracy : float
        Default 1e-3. Level of flux conservation during rescaling of kernel
    """
    def __init__(self, **kwargs):
        # sub_pixel_flag and flux_accuracy are taken care of in PSF base class
        super(FieldVaryingPSF, self).__init__(**kwargs)

        self.required_keys = ["filename"]
        utils.check_keys(self.meta, self.required_keys, action="error")

        self.meta["z_order"] = [61, 361]
        self._waveset, self.kernel_indexes = get_psf_wave_exts(self._file)
        self.current_ext = None
        self.current_data = None
        self._strehl_imagehdu = None

    def apply_to(self, fov):
        # .. todo: add in field rotation
        # accept "full", "dit", "none

        # check if there are any fov.fields to apply a psf to
        if len(fov.fields) > 0:
            if fov.hdu.data is None:
                fov.view(self.meta["sub_pixel_flag"])

            old_shape = fov.hdu.data.shape

            # get the kernels that cover this fov, and their respective masks
            # kernels and masks are returned by .get_kernel as a list of tuples
            canvas = None
            kernels_masks = self.get_kernel(fov)
            for kernel, mask in kernels_masks:

                # renormalise the kernel if needs be
                sum_kernel = np.sum(kernel)
                if abs(sum_kernel - 1) > self.meta["flux_accuracy"]:
                    kernel /= sum_kernel

                # image convolution
                new_image = convolve(fov.hdu.data, kernel, mode="same")
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
            fov.hdu.data = canvas

            # ..todo: careful with which dimensions mean what
            if "CRPIX1" in fov.hdu.header:
                fov.hdu.header["CRPIX1"] += (new_shape[0] - old_shape[0]) / 2
                fov.hdu.header["CRPIX2"] += (new_shape[1] - old_shape[1]) / 2

            if "CRPIX1D" in fov.hdu.header:
                fov.hdu.header["CRPIX1D"] += (new_shape[0] - old_shape[0]) / 2
                fov.hdu.header["CRPIX2D"] += (new_shape[1] - old_shape[1]) / 2

        return fov

    # def fov_grid(self, which="waveset"):
    #   This is taken care of by the PSF base class

    def get_kernel(self, fov):
        # 0. get file extension
        # 1. pull out strehl map for fov header
        # 2. get number of unique psfs
        # 3. pull out those psfs
        # 4. if more than one, make masks for the fov on the fov pixel scale
        # 5. make list of tuples with kernel and mask

        # find which file extension to use - keep a pointer in self.current_data
        fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        jj = nearest_index(fov_wave, self._waveset)
        ext = self.kernel_indexes[jj]
        if ext != self.current_ext:
            self.current_ext = ext
            self.current_data = self._file[ext].data

        # compare the fov and psf pixel scales
        kernel_pixel_scale = self._file[ext].header["CDELT1"]
        fov_pixel_scale = fov.hdu.header["CDELT1"]

        # get the spatial map of the kernel cube layers
        strl_hdu = self.strehl_imagehdu
        strl_cutout = get_strehl_cutout(fov.hdu.header, strl_hdu)

        # get the kernels and mask that fit inside the fov boundaries
        layer_ids = np.round(np.unique(strl_cutout.data)).astype(int)
        if len(layer_ids) > 1:
            kernels = [self.current_data[ii] for ii in layer_ids]
            # .. todo:: investigate. There's a .T in here that I don't like
            masks = [strl_cutout.data.T == ii for ii in layer_ids]
            self.kernel = [[krnl, msk] for krnl, msk in zip(kernels, masks)]
        else:
            self.kernel = [[self.current_data[layer_ids[0]], None]]

        # .. todo: re-scale kernel and masks to pixel_scale of FOV
        # .. todo: can this be put somewhere else to save on iterations?
        # .. todo: should the mask also be rescaled?
        # rescale the pixel scale of the kernel to match the fov images
        pix_ratio = fov_pixel_scale / kernel_pixel_scale
        if abs(pix_ratio - 1) > self.meta["flux_accuracy"]:
            for ii in range(len(self.kernel)):
                self.kernel[ii][0] = resize_array(self.kernel[ii][0], pix_ratio)

        return self.kernel

    @property
    def strehl_imagehdu(self):
        """ The HDU containing the positional info for kernel layers """
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


def resize_array(image, scale_factor, order=None):
    if order is None:
        order = rc.__currsys__["!SIM.computing.spline_order"]
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
