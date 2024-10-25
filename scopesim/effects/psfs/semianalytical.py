# -*- coding: utf-8 -*-
"""Currently only contains the AnisoCADO connection."""
from warnings import warn

import numpy as np
from scipy.interpolate import RectBivariateSpline
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

import anisocado as aniso

from .. import ter_curves_utils as tu
from ...base_classes import FieldOfViewBase
from ...utils import (figure_factory, figure_grid_factory, from_currsys,
                      quantify, check_keys)
from . import PSF


class SemiAnalyticalPSF(PSF):
    """Base class for semianalytical PSFs."""

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

        if self._kernel is not None:
            return self._kernel

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

        return self._kernel

    def remake_kernel(self, x):
        """
        Remake the kernel based on either a pixel_scale of FieldOfView.

        Parameters
        ----------
        x: float, FieldOfView
            [um] if float

        """
        warn("The 'remake_kernel' method was unused and thus deprecated and "
             "will be removed in a future release. If you are using this "
             "method, pleas let us know by creating an issue at: "
             "https://github.com/AstarVienna/ScopeSim/issues",
             DeprecationWarning, stacklevel=2)
        self._kernel = None
        return self.get_kernel(x)

    @property
    def wavelength(self):
        # FIXME: expensive property...
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
        nm_rms = nmrms_from_strehl_and_wavelength(strehl, wave, hdu)

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


class AnisocadoFieldVaryingPSF(AnisocadoConstPSF):

    required_keys = {"filename", "strehl", "wavelength", "r_max"}
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        params = {"psf_side_length": 256,
                  "grid_type": "radial",
                  "r_max": 30,
                  "anisocado_kwargs"}
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        self.psf_object = aniso.AnalyticalScaoPsf(pixelSize=pixel_scale,
                                                  N=self.meta["psf_side_length"],
                                                  wavelength=self.meta["wavelength"],
                                                  nmRms=self.nmRms,
                                                  **self.meta["anisocado_kwargs"])

        xs, ys = make_xys(r_max=self.meta["r_max"],
                          n=self.meta["psf_side_length"],
                          grid_type=self.meta["grid_type"])
        dec_idx_map = make_decimal_index_map(xs, ys, n=self.meta["psf_side_length"])
        self.psf_cube = generate_anisocado_psf_cube(xs, ys,
                                                    wave=self.meta["wavelength"],
                                                    pixel_scale=pixel_scale,
                                                    nmRms=self.nmRms,
                                                    psf_object=self.psf_object)


    def apply_to(self, obj, **kwargs):
        pass

    def get_kernel(self, fov):









def nmrms_from_strehl_and_wavelength(strehl: float,
                                     wavelength: float,
                                     strehl_hdu: fits.ImageHDU,
                                     plot: bool = False) -> float:
    """
    Return the wavefront error needed to make a PSF with desired strehl ratio.

    Parameters
    ----------
    strehl : float
        [0.001, 1] Desired strehl ratio. Values 1<sr<100 will be scale to <1
    wavelength : float
        [um]
    strehl_hdu : fits.ImageHDU
        2D map of strehl ratio as a function of wavelength [um] and residual
        wavefront error [nm RMS]
    plot : bool

    Returns
    -------
    nm : float
        [nm] residual wavefront error for generating an on-axis AnisoCADO PSF
        with the desired strehl ratio at a given wavelength

    """
    if 1. < strehl < 100.:
        strehl *= 0.01

    nm0 = strehl_hdu.header["CRVAL1"]
    dnm = strehl_hdu.header["CDELT1"]
    nm1 = strehl_hdu.header["NAXIS1"] * dnm + nm0
    nms = np.arange(nm0, nm1, dnm)

    w0 = strehl_hdu.header["CRVAL2"]
    dw = strehl_hdu.header["CDELT2"]
    w1 = strehl_hdu.header["NAXIS2"] * dw + w0
    ws = np.arange(w0, w1, dw)

    nms_spline = RectBivariateSpline(ws, nms, strehl_hdu.data, kx=1, ky=1)
    strehls = nms_spline(wavelength, nms)[0]

    if strehl > np.max(strehls):
        raise ValueError(f"Strehl ratio ({strehl}) is impossible at this "
                         f"wavelength ({wavelength}). Maximum Strehl possible "
                         f"is {np.max(strehls)}.")

    if strehls[0] < strehls[-1]:
        nm = np.interp(strehl, strehls, nms)
    else:
        nm = np.interp(strehl, strehls[::-1], nms[::-1])

    if plot:
        fig, ax = figure_factory()
        ax.plot(nms, strehls)
        ax.plot(nm, strehl, "ro")
        fig.show()

    return nm



def make_xys(r_max: float, n: int, grid_type: str="radial"):
    """
    Make a set of coordinates in the FIRST QUADRANT where PSFs can be generated

    .. note:: Only coordinates in the first quadrant are produced.
       Later functions will mirror these coordinates into the other 3 quadrants

    Parameters
    ----------
    r_max : float
        [arcsec] Side length of grid of points
    n : int
        Number of coords per side, e.g. `np.linspace(0, r_max, n)`
    grid_type : str, optional
        Default: "radial". Distribution of points in the grid.
        Options: ["square", "radial"]

    Returns
    -------
    xs, ys : float
        Two 1D arrays containing the x and y coordinates for the PSFs

    """

    assert grid_type in ["square", "radial"]
    if grid_type == "square":
        coords = [[x, y]
                  for x in np.linspace(0, r_max, n)
                  for y in np.linspace(0, r_max, n)]
    else:
        coords = [[r * np.cos(th), r * np.sin(th)]
                  for r in np.linspace(0, r_max, n)[1:]
                  for th in np.linspace(0, np.pi / 2, n)]
        coords = [[0, 0]] + coords

    coords = np.array(coords)
    xs, ys = coords[:, 0], coords[:, 1]

    return xs, ys


def make_decimal_index_map(xs, ys, n):
    """
    A PSF index map covering all 4 quadrants, decimals indicate PSF flipping

    The decimal index map covers the whole on-sky area where the PSFs are to
    be used. The assumption here is that the PSFs in quadrants 2,3,4 are the
    same as PSFs in the first quadrants, just flipped around either (or both) of
    the x and y axes. This assumption speeds up the calculation of the PSFs
    by a factor of 4

    The int(IMG) values give the index of the layer in the PSF cube for the
    relevant PSF for that area on sky. The decimal portion of the IMG values
    indicate the level of flipping needed for the PSF kernel image

    Quadrant fractions are defined as:
    - x.00 = Quadrant 1 : (+, +) -> leave as is
    - x.25 = Quadrant 2 : (-, +) ->  np.fliplr
    - x.50 = Quadrant 4 : (+, -) -> np.flipud
    - x.75 = Quadrant 3 : (-, -) -> np.fliplr and np.flipud

    Parameters
    ----------
    xs, ys : array-like
        [arcsec] Coordinated where PSFs should be generated
    n : int
        Side length of the PSF index map

    Returns
    -------
    hdu : fits.ImageHDU
        The ImageHDU containing the decimal index map


    """
    assert len(xs) == len(ys)

    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()

    x_grid = np.linspace(x_min, x_max, n)
    y_grid = np.linspace(y_min, y_max, n)
    X, Y = np.meshgrid(x_grid, y_grid)

    # Make a meshgrid of the distances of each pixel from the coordinate of each PSF
    X_dists = X[None, :, :] - xs[:, None, None]
    Y_dists = Y[None, :, :] - ys[:, None, None]
    R_dists = np.sqrt(X_dists ** 2 + Y_dists ** 2)
    # Sort the cube of distances and take the closest index as the PSF that is relevant for each pixel
    I = np.argsort(R_dists, axis=0)[0]

    # Add the other 3 quadrants as h- or v-flipped versions of the original index map
    I2 = np.hstack([np.fliplr(I) + 0.25, I])
    I3 = np.vstack([np.flipud(I2) + 0.5, I2])

    # Create a FITS ImageHDU with the I3 map and header WCS
    w = WCS()
    w.cdelt = [x_grid[1] - x_grid[0], y_grid[1] - y_grid[0]]
    w.cunit = [u.arcsec, u.arcsec]
    w.ctype = ["LINEAR", "LINEAR"]
    w.crval = [0, 0]
    w.crpix = [n, n]

    hdu = fits.ImageHDU(header=w.to_header(), data=I3)

    return hdu


def generate_anisocado_psf_cube(dxs, dys, wave=2.15,
                                pixel_scale=0.004, n=256, nmRms=300,
                                psf_object=None, **anisocado_kwargs):
    """
    Creates a cube of off-axis SCAO PSF kernels using AnisoCADO

    Parameters
    ----------
    dxs, dys : array-like
        [arcsec] Positions relative to the optical axis for PSF kernels
    wave : float, optional
        [um] Wavelength of PSF kernel. Default 2.15um
    pixel_scale : float, optional
        [arcsec] Default 0.004" for MICADO
    n : int, optional
        [pixel] Default 256. Side length of PSF kernel image
    nmRms : float, optional
        [nm] Residual wavefront error. Generally taken from a Strehl-nmRms map
    anisocado_kwargs : dict
        Additional kwargs to be passed directly to anisocado.AnalyticalScaoPsf

    Returns
    -------

    """
    if psf_object is None:
        psf_object = aniso.AnalyticalScaoPsf(pixelSize=pixel_scale,
                                             N=n,
                                             wavelength=wave,
                                             nmRms=nmRms,
                                             **anisocado_kwargs)
    # ..todo : check orientation of PSF in real MICADO images re ' .T'
    psf_cube = np.array([psf_object.shift_off_axis(dx, dy).T
                         for dx, dy in zip(dxs, dys)])

    # round the edges
    def round_edges(kernel):
        """
        Sets all pixels below a threshold value to 0 to 'round' out the kernel
        """
        x, y = np.array(kernel.shape) // 2
        z = np.min([kernel[x, 0], kernel[x, -1], kernel[0, y], kernel[-1, y]])
        kernel[kernel < z] = 0.
        return kernel

    psf_cube = np.array([round_edges(kernel) for kernel in psf_cube])

    # Create a FITS ImageHDU with the I3 map and header WCS
    w = WCS()
    w.cdelt = [pixel_scale, pixel_scale]
    w.cunit = [u.arcsec, u.arcsec]
    w.ctype = ["LINEAR", "LINEAR"]
    w.crval = [0, 0]
    w.crpix = [n / 2, n / 2]

    hdu = fits.ImageHDU(header=w.to_header(), data=psf_cube)
    hdu.header["WAVE0"] = wave
    hdu.header["NM_RMS"] = nmRms

    return psf_cube


def get_psf_from_decimal_index(psf_cube, dec_idx):
    """
    Flips the PSF kernel according to the quadrant where it should be applied

    Quadrant fractions are defined as:
    - x.00 = Quadrant 1 : (+, +) -> leave as is
    - x.25 = Quadrant 2 : (-, +) ->  np.fliplr
    - x.50 = Quadrant 4 : (+, -) -> np.flipud
    - x.75 = Quadrant 3 : (-, -) -> np.flipud

    Parameters
    ----------
    psf_cube : np.ndarray
        3D array of PSF kernels
    dec_idx : float
        Decimal index = [PSF index].[quadrant fraction]

    Returns
    -------
    psf_im : np.ndarray
        2D array of PSF kernel for the given decimal index

    """
    assert dec_idx < len(psf_cube)

    psf_im = psf_cube[int(dec_idx)]
    frac_idx = dec_idx - int(dec_idx)

    if frac_idx >= 0.5:
        psf_im = np.flipud(psf_im)
        frac_idx -= 0.5
    if frac_idx > 0.:
        psf_im = np.fliplr(psf_im)

    return psf_im
