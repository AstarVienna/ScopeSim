# -*- coding: utf-8 -*-
"""."""

import warnings

import numpy as np
from astropy import units as u
from astropy.convolution import Gaussian2DKernel

from ...base_classes import ImagePlaneBase, FieldOfViewBase
from ...utils import from_currsys, quantify, figure_factory, check_keys
from . import PSF, PoorMansFOV
from . import psf_utils as pu


class AnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [41, 641]
        self.convolution_classes = FieldOfViewBase


class Vibration(AnalyticalPSF):
    """Creates a wavelength independent kernel image."""

    required_keys = {"fwhm", "pixel_scale"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [244, 744]
        self.meta["width_n_fwhms"] = 4
        self.convolution_classes = ImagePlaneBase

        check_keys(self.meta, self.required_keys, action="error")
        self.kernel = None

    def get_kernel(self, obj):
        if self.kernel is None:
            from_currsys(self.meta, self.cmds)
            fwhm_pix = self.meta["fwhm"] / self.meta["pixel_scale"]
            sigma = fwhm_pix / 2.35
            width = max(1, int(fwhm_pix * self.meta["width_n_fwhms"]))
            self.kernel = Gaussian2DKernel(sigma, x_size=width, y_size=width,
                                           mode="center").array
            self.kernel /= np.sum(self.kernel)

        return self.kernel.astype(float)


class NonCommonPathAberration(AnalyticalPSF):
    """
    TBA.

    Needed: pixel_scale
    Accepted: kernel_width, strehl_drift
    """

    required_keys = {"pixel_scale"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [241, 641]
        self.meta["kernel_width"] = None
        self.meta["strehl_drift"] = 0.02
        self.meta["wave_min"] = "!SIM.spectral.wave_min"
        self.meta["wave_max"] = "!SIM.spectral.wave_max"

        self._total_wfe = None

        self.valid_waverange = [0.1 * u.um, 0.2 * u.um]

        self.convolution_classes = FieldOfViewBase
        check_keys(self.meta, self.required_keys, action="error")

    def fov_grid(self, which="waveset", **kwargs):
        """See parent docstring."""
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release.", DeprecationWarning, stacklevel=2)
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = from_currsys(self.meta, self.cmds)

            min_sr = pu.wfe2strehl(self.total_wfe, self.meta["wave_min"])
            max_sr = pu.wfe2strehl(self.total_wfe, self.meta["wave_max"])

            srs = np.arange(min_sr, max_sr, self.meta["strehl_drift"])
            waves = 6.2831853 * self.total_wfe * (-np.log(srs))**-0.5
            waves = quantify(waves, u.um).value
            waves = (list(waves) + [self.meta["wave_max"]]) * u.um
        else:
            waves = [] * u.um

        return waves

    def get_kernel(self, obj):
        waves = obj.meta["wave_min"], obj.meta["wave_max"]

        old_waves = self.valid_waverange
        wave_mid_old = 0.5 * (old_waves[0] + old_waves[1])
        wave_mid_new = 0.5 * (waves[0] + waves[1])
        strehl_old = pu.wfe2strehl(wfe=self.total_wfe, wave=wave_mid_old)
        strehl_new = pu.wfe2strehl(wfe=self.total_wfe, wave=wave_mid_new)

        if np.abs(1 - strehl_old / strehl_new) > self.meta["strehl_drift"]:
            self.valid_waverange = waves
            self.kernel = pu.wfe2gauss(wfe=self.total_wfe, wave=wave_mid_new,
                                       width=self.meta["kernel_width"])
            self.kernel /= np.sum(self.kernel)

        return self.kernel

    @property
    def total_wfe(self):
        if self._total_wfe is None:
            if self.table is not None:
                self._total_wfe = pu.get_total_wfe_from_table(self.table)
            else:
                self._total_wfe = 0
        return self._total_wfe

    def plot(self):
        fig, axes = figure_factory()

        wave_min, wave_max = from_currsys([self.meta["wave_min"],
                                           self.meta["wave_max"]], self.cmds)
        waves = np.linspace(wave_min, wave_max, 1001) * u.um
        wfe = self.total_wfe
        strehl = pu.wfe2strehl(wfe=wfe, wave=waves)

        axes.plot(waves, strehl)
        axes.set_xlabel(f"Wavelength [{waves.unit}]")
        axes.set_ylabel(f"Strehl Ratio \n[Total WFE = {wfe}]")

        return fig


class SeeingPSF(AnalyticalPSF):
    """
    Currently only returns gaussian kernel with a ``fwhm`` [arcsec].

    Parameters
    ----------
    fwhm : flaot
        [arcsec]

    """

    def __init__(self, fwhm=1.5, **kwargs):
        super().__init__(**kwargs)

        self.meta["fwhm"] = fwhm
        self.meta["z_order"] = [242, 642]

    def get_kernel(self, fov):
        # called by .apply_to() from the base PSF class

        pixel_scale = fov.header["CDELT1"] * u.deg.to(u.arcsec)
        pixel_scale = quantify(pixel_scale, u.arcsec)

        # add in the conversion to fwhm from seeing and wavelength here
        fwhm = from_currsys(self.meta["fwhm"], self.cmds) * u.arcsec / pixel_scale

        sigma = fwhm.value / 2.35
        kernel = Gaussian2DKernel(sigma, mode="center").array
        kernel /= np.sum(kernel)

        return kernel

    def plot(self):
        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        spec_dict = from_currsys("!SIM.spectral", self.cmds)
        return super().plot(PoorMansFOV(pixel_scale, spec_dict))


class GaussianDiffractionPSF(AnalyticalPSF):
    def __init__(self, diameter, **kwargs):
        super().__init__(**kwargs)
        self.meta["diameter"] = diameter
        self.meta["z_order"] = [242, 642]

    def fov_grid(self, which="waveset", **kwargs):
        """See parent docstring."""
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release.", DeprecationWarning, stacklevel=2)
        wavelengths = []
        if which == "waveset" and \
                "waverange" in kwargs and \
                "pixel_scale" in kwargs:
            waverange = quantify(kwargs["waverange"], u.um)
            diameter = quantify(self.meta["diameter"], u.m).to(u.um)
            fwhm = 1.22 * (waverange / diameter).value  # in rad

            pixel_scale = quantify(kwargs["pixel_scale"], u.deg)
            pixel_scale = pixel_scale.to(u.rad).value
            fwhm_range = np.arange(fwhm[0], fwhm[1], pixel_scale)
            wavelengths = list(fwhm_range / 1.22 * diameter.to(u.m))

        # TODO: check that this is actually correct
        return wavelengths

    def update(self, **kwargs):
        if "diameter" in kwargs:
            self.meta["diameter"] = kwargs["diameter"]

    def get_kernel(self, fov):
        # called by .apply_to() from the base PSF class

        pixel_scale = fov.header["CDELT1"] * u.deg.to(u.arcsec)
        pixel_scale = quantify(pixel_scale, u.arcsec)

        wave = 0.5 * (fov.meta["wave_max"] + fov.meta["wave_min"])

        wave = quantify(wave, u.um)
        diameter = quantify(self.meta["diameter"], u.m).to(u.um)
        fwhm = 1.22 * (wave / diameter) * u.rad.to(u.arcsec) / pixel_scale

        sigma = fwhm.value / 2.35
        kernel = Gaussian2DKernel(sigma, mode="center").array
        kernel /= np.sum(kernel)

        return kernel

    def plot(self):
        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        spec_dict = from_currsys("!SIM.spectral", self.cmds)
        return super().plot(PoorMansFOV(pixel_scale, spec_dict))
