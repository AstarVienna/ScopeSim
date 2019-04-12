import numpy as np
from astropy import units as u
from astropy.convolution import Gaussian2DKernel

from ... import utils
from .psfs import AnalyticalPSF


class GaussianDiffractionPSF(AnalyticalPSF):
    def __init__(self, diameter, **kwargs):
        super(GaussianDiffractionPSF, self).__init__(**kwargs)
        self.meta["diameter"] = diameter
        self.meta["z_order"] = [0, 300]

    def fov_grid(self, header=None, waverange=None, **kwargs):
        waverange = utils.quantify(waverange, u.um)
        diameter = utils.quantify(self.meta["diameter"], u.m).to(u.um)
        fwhm = 1.22 * (waverange / diameter).value  # in rad

        pixel_scale = utils.quantify(header["CDELT1"], u.deg).to(u.rad).value
        fwhm_range = np.arange(fwhm[0], fwhm[1], pixel_scale)
        wavelengths = fwhm_range / 1.22 * diameter.to(u.m)

        return {"coords": None, "wavelengths": wavelengths}

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
