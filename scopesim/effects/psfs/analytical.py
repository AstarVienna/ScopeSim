# -*- coding: utf-8 -*-
"""Contains simple Vibration, NCPA, Seeing and Diffraction PSFs."""

from typing import ClassVar

import numpy as np
from astropy import units as u
from astropy.convolution import Gaussian2DKernel

from ...base_classes import ImagePlaneBase, FieldOfViewBase
from ...utils import (from_currsys, quantify, quantity_from_table,
                      figure_factory, check_keys)
from . import PSF, PoorMansFOV


class AnalyticalPSF(PSF):
    """Base class for analytical PSFs."""

    z_order: ClassVar[tuple[int, ...]] = (41, 641)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.convolution_classes = FieldOfViewBase


class Vibration(AnalyticalPSF):
    """Creates a wavelength independent kernel image."""

    required_keys = {"fwhm", "pixel_scale"}
    z_order: ClassVar[tuple[int, ...]] = (244, 744)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["width_n_fwhms"] = 4
        self.convolution_classes = ImagePlaneBase

        check_keys(self.meta, self.required_keys, action="error")
        self.kernel = None

    def get_kernel(self, obj):
        if self.kernel is not None:
            return self.kernel

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
    z_order: ClassVar[tuple[int, ...]] = (241, 641)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["kernel_width"] = None
        self.meta["strehl_drift"] = 0.02
        self.meta["wave_min"] = "!SIM.spectral.wave_min"
        self.meta["wave_max"] = "!SIM.spectral.wave_max"

        self._total_wfe = None

        self.valid_waverange = [0.1 * u.um, 0.2 * u.um]

        self.convolution_classes = FieldOfViewBase
        check_keys(self.meta, self.required_keys, action="error")

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
            self.kernel /= np.sum(self.kernel)

        return self.kernel

    def _get_total_wfe_from_table(self):
        wfes = quantity_from_table("wfe_rms", self.table, "um")
        n_surfs = self.table["n_surfaces"]
        return np.sum(n_surfs * wfes**2)**0.5

    @property
    def total_wfe(self):
        if self._total_wfe is not None:
            return self._total_wfe

        if self.table is not None:
            self._total_wfe = self._get_total_wfe_from_table()
        else:
            self._total_wfe = 0

        return self._total_wfe

    def plot(self):
        fig, axes = figure_factory()

        wave_min, wave_max = from_currsys([self.meta["wave_min"],
                                           self.meta["wave_max"]], self.cmds)
        waves = np.linspace(wave_min, wave_max, 1001) * u.um
        wfe = self.total_wfe
        strehl = wfe2strehl(wfe=wfe, wave=waves)

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

    z_order: ClassVar[tuple[int, ...]] = (242, 642)

    def __init__(self, fwhm=1.5, **kwargs):
        super().__init__(**kwargs)

        self.meta["fwhm"] = fwhm

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
    z_order: ClassVar[tuple[int, ...]] = (242, 642)

    def __init__(self, diameter, **kwargs):
        super().__init__(**kwargs)
        self.meta["diameter"] = diameter

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


def wfe2gauss(wfe, wave, width=None):
    strehl = wfe2strehl(wfe, wave)
    sigma = _strehl2sigma(strehl)
    if width is None:
        width = int(np.ceil(8 * sigma))
        width += (width + 1) % 2
    gauss = _sigma2gauss(sigma, x_size=width, y_size=width)

    return gauss


def wfe2strehl(wfe, wave):
    wave = quantify(wave, u.um)
    wfe = quantify(wfe, u.um)
    x = 2 * 3.1415926526 * wfe / wave
    strehl = np.exp(-x**2)
    return strehl


def _strehl2sigma(strehl):
    amplitudes = [0.00465, 0.00480, 0.00506, 0.00553, 0.00637, 0.00793,
                  0.01092, 0.01669, 0.02736, 0.04584, 0.07656, 0.12639,
                  0.20474, 0.32156, 0.48097, 0.66895, 0.84376, 0.95514,
                  0.99437, 0.99982, 0.99999]
    sigmas = [19.9526, 15.3108, 11.7489, 9.01571, 6.91830, 5.30884, 4.07380,
              3.12607, 2.39883, 1.84077, 1.41253, 1.08392, 0.83176, 0.63826,
              0.48977, 0.37583, 0.28840, 0.22130, 0.16982, 0.13031, 0.1]
    sigma = np.interp(strehl, amplitudes, sigmas)
    return sigma


def _sigma2gauss(sigma, x_size=15, y_size=15):
    kernel = Gaussian2DKernel(sigma, x_size=x_size, y_size=y_size,
                              mode="oversample").array
    kernel /= np.sum(kernel)
    return kernel
