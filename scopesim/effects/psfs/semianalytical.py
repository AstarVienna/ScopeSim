# -*- coding: utf-8 -*-
"""."""

import numpy as np
from astropy import units as u

import anisocado as aniso

from .. import ter_curves_utils as tu
from ...base_classes import FieldOfViewBase
from ...utils import figure_grid_factory, from_currsys, quantify, check_keys
from . import PSF
from . import psf_utils as pu


class SemiAnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [42]
        self.convolution_classes = FieldOfViewBase
        # self.convolution_classes = ImagePlaneBase


class AnisocadoConstPSF(SemiAnalyticalPSF):
    """
    Makes a SCAO on-axis PSF with a desired Strehl ratio at a given wavelength.

    To make the PSFs a map connecting Strehl, Wavelength, and residual
    wavefront error is required.

    Parameters
    ----------
    filename : str
        Path to Strehl map with axes (x, y) = (wavelength, wavefront error).
    strehl : float
        Desired Strehl ratio. Either percentage [1, 100] or fractional
        [1e-3, 1].
    wavelength : float
        [um] The given strehl is valid for this wavelength.
    psf_side_length : int
        [pixel] Default is 512. Side length of the kernel images.
    offset : tuple
        [arcsec] SCAO guide star offset from centre (dx, dy).
    rounded_edges : bool
        Default is True. Sets all halo values below a threshold to zero.
        The threshold is determined from the max values of the edge rows of the
        kernel image.

    Other Parameters
    ----------------
    convolve_mode : str
        ["same", "full"] convolution keywords from scipy.signal.convolve

    Examples
    --------
    Add an AnisocadoConstPSF with code::

        from scopesim.effects import AnisocadoConstPSF
        psf = AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                                strehl=0.5,
                                wavelength=2.15,
                                convolve_mode="same",
                                psf_side_length=512)

    Add an AnisocadoConstPSF to a yaml file::

        effects:
        -   name: Ks_Stehl_40_PSF
            description: A 40% Strehl PSF over the field of view
            class: AnisocadoConstPSF
            kwargs:
                filename: "test_AnisoCADO_rms_map.fits"
                strehl: 0.5
                wavelength: 2.15
                convolve_mode: full
                psf_side_length: 512

    """

    required_keys = {"filename", "strehl", "wavelength"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "z_order": [42, 652],
            "psf_side_length": 512,
            "offset": (0, 0),
            "rounded_edges": True,
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")
        self.nmRms      # check to see if it throws an error

        self._psf_object = None
        self._kernel = None

    def get_kernel(self, fov):
        # called by .apply_to() from the base PSF class

        if self._kernel is None:
            if isinstance(fov, FieldOfViewBase):
                pixel_scale = fov.header["CDELT1"] * u.deg.to(u.arcsec)
            elif isinstance(fov, float):
                pixel_scale = fov

            n = self.meta["psf_side_length"]
            wave = self.wavelength
            self._psf_object = aniso.AnalyticalScaoPsf(pixelSize=pixel_scale,
                                                       N=n, wavelength=wave,
                                                       nmRms=self.nmRms)
            if np.any(self.meta["offset"]):
                self._psf_object.shift_off_axis(self.meta["offset"][0],
                                                self.meta["offset"][1])

            self._kernel = self._psf_object.psf_latest
            self._kernel /= np.sum(self._kernel)
            if self.meta["rounded_edges"]:
                self._kernel = pu.round_kernel_edges(self._kernel)

        return self._kernel

    def remake_kernel(self, x):
        """
        Remake the kernel based on either a pixel_scale of FieldOfView.

        Parameters
        ----------
        x: float, FieldOfView
            [um] if float

        """
        self._kernel = None
        return self.get_kernel(x)

    @property
    def wavelength(self):
        wave = from_currsys(self.meta["wavelength"], self.cmds)
        if isinstance(wave, str) and wave in tu.FILTER_DEFAULTS:
            filter_name = from_currsys(wave, cmds=self.cmds)
            wave = tu.get_filter_effective_wavelength(filter_name)
        wave = quantify(wave, u.um).value

        return wave

    @property
    def strehl_ratio(self):
        strehl = None
        if self._psf_object is not None:
            strehl = self._psf_object.strehl_ratio

        return strehl

    @property
    def nmRms(self):
        strehl = from_currsys(self.meta["strehl"], self.cmds)
        wave = self.wavelength
        hdu = self._file[0]
        nm_rms = pu.nmrms_from_strehl_and_wavelength(strehl, wave, hdu)

        return nm_rms

    def plot(self, obj=None, **kwargs):
        fig, gs = figure_grid_factory(
            2, 2, height_ratios=(3, 2),
            left=0.3, right=0.7, bottom=0.15, top=0.85,
            wspace=0.05, hspace=0.05)
        # or no height_ratios and bottom=0.1, top=0.9

        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        kernel = self.get_kernel(pixel_scale)

        ax = fig.add_subplot(gs[0, 0])
        im = kernel
        r_sky = pixel_scale * im.shape[0]
        ax.imshow(im, norm="log", origin="lower",
                  extent=[-r_sky, r_sky, -r_sky, r_sky], **kwargs)
        ax.set_aspect("equal")
        ax.set_xlabel("[arcsec]")
        ax.set_ylabel("[arcsec]")
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")

        ax = fig.add_subplot(gs[0, 1])
        x = kernel.shape[1] // 2
        y = kernel.shape[0] // 2
        r = 16
        im = kernel[y-r:y+r, x-r:x+r]
        r_sky = pixel_scale * im.shape[0]
        ax.imshow(im, norm="log", origin="lower",
                  extent=[-r_sky, r_sky, -r_sky, r_sky], **kwargs)
        ax.set_aspect("equal")
        ax.set_xlabel("[arcsec]")
        ax.set_ylabel("[arcsec]")
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("right")

        ax = fig.add_subplot(gs[1, :])
        hdr = self._file[0].header
        data = self._file[0].data
        wfes = np.arange(hdr["NAXIS1"]) * hdr["CDELT1"] + hdr["CRVAL1"]
        waves = np.arange(hdr["NAXIS2"]) * hdr["CDELT2"] + hdr["CRVAL2"]

        # TODO: Get unit dynamically? Then again, it's hardcoded elsewhere in
        #       this module...
        unit_str = u.Unit("um").to_string("latex")
        for strehl, wav in reversed(list(zip(data, waves))):
            ax.plot(wfes, strehl, label=f"{wav:.3f} {unit_str}")

        ax.set_xlabel(f"RMS Wavefront Error [{unit_str}]")
        ax.set_ylabel("Strehl Ratio")
        ax.legend()
        fig.align_labels()
        return fig
