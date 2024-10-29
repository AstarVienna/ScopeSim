# -*- coding: utf-8 -*-
"""
Utility classes and functions for ``SpectralTraceList``.

This module contains
   - the definition of the `SpectralTrace` class. The visible effect should
     always be a `SpectralTraceList`, even if that contains only one
     `SpectralTrace`.
   - the definition of the `XiLamImage` class
   - utility functions for use with spectral traces
"""

import warnings

import numpy as np

from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d

from astropy.table import Table, vstack
from astropy.modeling import fitting
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.modeling.models import Polynomial2D

from ..utils import (power_vector, quantify, from_currsys, close_loop,
                     figure_factory, get_logger)


logger = get_logger(__name__)


class SpectralTrace:
    """Definition of one spectral trace.

    A SpectralTrace describes the mapping of spectral slit coordinates
    to the focal plane. The class reads an order layout and fits several
    functions to describe the geometry of the trace.

    Slit coordinates are:
    - xi : spatial position along the slit [arcsec]
    - lam : Wavelength [um]
    Focal plane coordinates are:
    - x, y : [mm]
    """

    _class_params = {
        "x_colname": "x",
        "y_colname": "y",
        "s_colname": "s",
        "wave_colname": "wavelength",
        "dwave": 0.002,
        "aperture_id": 0,
        "image_plane_id": 0,
        "extension_id": 2,
        "spline_order": 4,
        "pixel_size": None,
        "description": "<no description>",
    }

    def __init__(self, trace_tbl, cmds=None, **kwargs):
        # Within scopesim, the actual parameter values are
        # passed as kwargs from the SpectralTraceList.
        # The values need to be here for stand-alone use.
        self.meta = {}
        self.meta.update(self._class_params)
        self.meta.update(kwargs)
        self.cmds = cmds

        if isinstance(trace_tbl, (fits.BinTableHDU, fits.TableHDU)):
            self.table = Table.read(trace_tbl)
            self.meta["trace_id"] = trace_tbl.header.get("EXTNAME",
                                                         "<unknown trace id>")
            self.dispersion_axis = trace_tbl.header.get("DISPDIR", "unknown")
        elif isinstance(trace_tbl, Table):
            self.table = trace_tbl
            self.dispersion_axis = "unknown"
        else:
            raise ValueError("trace_tbl must be one of (fits.BinTableHDU, "
                             "fits.TableHDU, astropy.Table) but is "
                             f"{type(trace_tbl)}")
        self.compute_interpolation_functions()

        # Declaration of other attributes
        self._xilamimg = None
        self.dlam_per_pix = None

    @property
    def trace_id(self):
        """Return the name of the trace."""
        return self.meta["trace_id"]

    def fov_grid(self):
        """
        Provide information on the source space volume required by the effect.

        Returns
        -------
        A dictionary with entries `wave_min` and `wave_max`.
        Spatial limits are determined by the `ApertureMask` effect
        and are not returned here.
        """
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release.", DeprecationWarning, stacklevel=2)
        aperture_id = self.meta["aperture_id"]
        lam_arr = self.table[self.meta["wave_colname"]]

        wave_max = np.max(lam_arr)
        wave_min = np.min(lam_arr)

        return {"wave_min": wave_min, "wave_max": wave_max,
                "trace_id": self.trace_id, "aperture_id": aperture_id}

    def compute_interpolation_functions(self):
        """
        Compute various interpolation functions between slit and focal plane.

        Focal plane coordinates are `x` and `y`, in mm. Slit coordinates are
        `xi` (spatial coordinate along the slit, in arcsec) and `lam`
        (wavelength, in um).
        """
        x_arr = self.table[self.meta["x_colname"]]
        y_arr = self.table[self.meta["y_colname"]]
        xi_arr = self.table[self.meta["s_colname"]]
        lam_arr = self.table[self.meta["wave_colname"]]

        self.wave_min = quantify(np.min(lam_arr), u.um).value
        self.wave_max = quantify(np.max(lam_arr), u.um).value

        self.xy2xi = Transform2D.fit(x_arr, y_arr, xi_arr)
        self.xy2lam = Transform2D.fit(x_arr, y_arr, lam_arr)
        self.xilam2x = Transform2D.fit(xi_arr, lam_arr, x_arr)
        self.xilam2y = Transform2D.fit(xi_arr, lam_arr, y_arr)
        self._xiy2x = Transform2D.fit(xi_arr, y_arr, x_arr)
        self._xiy2lam = Transform2D.fit(xi_arr, y_arr, lam_arr)

        if self.dispersion_axis == "unknown":
            dlam_dx, dlam_dy = self.xy2lam.gradient()
            wave_mid = 0.5 * (self.wave_min + self.wave_max)
            xi_mid = np.mean(xi_arr)
            x_mid = self.xilam2x(xi_mid, wave_mid)
            y_mid = self.xilam2y(xi_mid, wave_mid)
            if dlam_dx(x_mid, y_mid) > dlam_dy(x_mid, y_mid):
                self.dispersion_axis = "x"
            else:
                self.dispersion_axis = "y"
            logger.info(
                "Dispersion axis determined to be %s", self.dispersion_axis)

    def map_spectra_to_focal_plane(self, fov):
        """
        Apply the spectral trace mapping to a spectral cube.

        The cube is contained in a FieldOfView object, which also has
        world coordinate systems for the Source (sky coordinates and
        wavelengths) and for the focal plane.
        The method returns a section of the fov image along with info on
        where this image lies in the focal plane.
        """
        logger.info("Mapping %s", fov.trace_id)
        # Initialise the image based on the footprint of the spectral
        # trace and the focal plane WCS
        wave_min = fov.meta["wave_min"].value       # [um]
        wave_max = fov.meta["wave_max"].value       # [um]
        xi_min = fov.meta["xi_min"].value           # [arcsec]
        xi_max = fov.meta["xi_max"].value           # [arcsec]
        xlim_mm, ylim_mm = self.footprint(wave_min=wave_min, wave_max=wave_max,
                                          xi_min=xi_min, xi_max=xi_max)

        if xlim_mm is None:
            logger.warning("xlim_mm is None")
            return None

        fov_header = fov.header
        det_header = fov.detector_header

        # WCSD from the FieldOfView - this is the full detector plane
        pixsize = fov_header["CDELT1D"] * u.Unit(fov_header["CUNIT1D"])
        pixsize = pixsize.to(u.mm).value
        pixscale = fov_header["CDELT1"] * u.Unit(fov_header["CUNIT1"])
        pixscale = pixscale.to(u.arcsec).value

        fpa_wcsd = WCS(det_header, key="D")
        naxis1d, naxis2d = det_header["NAXIS1"], det_header["NAXIS2"]
        xlim_px, ylim_px = fpa_wcsd.all_world2pix(xlim_mm, ylim_mm, 0)
        xmin = np.floor(xlim_px.min()).astype(int)
        xmax = np.ceil(xlim_px.max()).astype(int)
        ymin = np.floor(ylim_px.min()).astype(int)
        ymax = np.ceil(ylim_px.max()).astype(int)

        # Check if spectral trace footprint is outside FoV
        if xmax < 0 or xmin > naxis1d or ymax < 0 or ymin > naxis2d:
            logger.info(
                "Spectral trace %s: footprint is outside FoV", fov.trace_id)
            return None

        # Only work on parts within the FoV
        xmin = max(xmin, 0)
        xmax = min(xmax, naxis1d)
        ymin = max(ymin, 0)
        ymax = min(ymax, naxis2d)

        # Create header for the subimage - I think this only needs the DET one,
        # but we'll do both. The WCSs are initialised from the full fpa WCS and
        # then shifted accordingly.
        det_wcs = WCS(det_header, key="D")
        det_wcs.wcs.crpix -= np.array([xmin, ymin])

        sub_naxis1 = xmax - xmin
        sub_naxis2 = ymax - ymin

        # initialise the subimage
        image = np.zeros((sub_naxis2, sub_naxis1), dtype=np.float32)

        # Adjust the limits of the subimage in millimeters in the focal plane
        # This takes the adjustment to integer pixels into account
        xmin_mm, ymin_mm = fpa_wcsd.all_pix2world(xmin, ymin, 0)
        xmax_mm, ymax_mm = fpa_wcsd.all_pix2world(xmax, ymax, 0)

        self._set_dispersion(wave_min, wave_max, pixsize=pixsize)
        try:
            xilam = XiLamImage(fov, self.dlam_per_pix)
            self._xilamimg = xilam   # ..todo: remove or make available with a debug flag?
        except ValueError:
            logger.warning(" ---> %s gave ValueError", self.trace_id)

        npix_xi, npix_lam = xilam.npix_xi, xilam.npix_lam
        xilam_wcs = xilam.wcs

        # focal-plane coordinate images
        ximg_fpa, yimg_fpa = np.meshgrid(np.linspace(xmin_mm, xmax_mm,
                                                     sub_naxis1,
                                                     dtype=np.float32),
                                         np.linspace(ymin_mm, ymax_mm,
                                                     sub_naxis2,
                                                     dtype=np.float32))

        # Image mapping (xi, lambda) on the focal plane
        xi_fpa = self.xy2xi(ximg_fpa, yimg_fpa).astype(np.float32)
        lam_fpa = self.xy2lam(ximg_fpa, yimg_fpa).astype(np.float32)

        # mask everything outside the wavelength range
        mask = (xi_fpa >= xi_min) & (xi_fpa <= xi_max)
        xi_fpa *= mask
        lam_fpa *= mask

        # mask everything outside the wavelength range
        # mask = (lam_fpa >= wave_min) & (lam_fpa <= wave_max)
        # xi_fpa *= mask
        # lam_fpa *= mask

        # Convert to pixel images
        # These are the pixel coordinates in the image corresponding to
        # xi, lambda
        # It is much quicker to do the linear transformation by hand
        # than to use the astropy.wcs functions for conversion.
        i_img = ((lam_fpa - xilam_wcs.wcs.crval[0])
                 / xilam_wcs.wcs.cdelt[0]).astype(int)
        j_img = ((xi_fpa - xilam_wcs.wcs.crval[1])
                 / xilam_wcs.wcs.cdelt[1]).astype(int)

        # truncate images to remove pixel coordinates outside the image
        ijmask = ((i_img >= 0) * (i_img < npix_lam)
                  * (j_img >= 0) * (j_img < npix_xi))

        # do the actual interpolation
        # image is in [ph/s/um/arcsec]
        image = xilam.interp(xi_fpa, lam_fpa, grid=False) * ijmask

        # Scale to ph / s / pixel
        dlam_by_dx, dlam_by_dy = self.xy2lam.gradient()
        dlam_per_pix = pixsize * np.sqrt(dlam_by_dx(ximg_fpa, yimg_fpa)**2 +
                                         dlam_by_dy(ximg_fpa, yimg_fpa)**2)
        image *= pixscale * dlam_per_pix        # [arcsec/pix] * [um/pix]

        # img_header = sub_wcs.to_header()
        # img_header.update(det_wcs.to_header())
        img_header = det_wcs.to_header()
        img_header["XMIN"] = xmin
        img_header["XMAX"] = xmax
        img_header["YMIN"] = ymin
        img_header["YMAX"] = ymax

        if np.any(image < 0):
            logger.warning("map_spectra_to_focal_plane: %d negative pixels",
                           np.sum(image < 0))

        image_hdu = fits.ImageHDU(header=img_header, data=image)
        return image_hdu

    def rectify(self, hdulist, interps=None, wcs=None, **kwargs):
        """Create 2D spectrum for a trace.

        Parameters
        ----------
        hdulist : HDUList
           The result of scopesim readout
        interps : list of interpolation functions
           If provided, there must be one for each image extension in `hdulist`.
           The functions go from pixels to the images and can be created with,
           e.g., RectBivariateSpline.
        wcs : The WCS describing the rectified XiLamImage. This can be created
           in a simple way from the fov included in the `OpticalTrain` used in
           the simulation run producing `hdulist`.

        The WCS can also be set up via the following keywords:

        bin_width : float [um]
           The spectral bin width. This is best computed automatically from the
           spectral dispersion of the trace.
        wave_min, wave_max : float [um]
           Limits of the wavelength range to extract. The default is the
           the full range on which the `SpectralTrace` is defined. This may
           extend significantly beyond the filter window.
        xi_min, xi_max : float [arcsec]
           Spatial limits of the slit on the sky. This should be taken from
           the header of the hdulist, but this is not yet provided by scopesim
        """
        logger.info("Rectifying %s", self.trace_id)

        wave_min = kwargs.get("wave_min",
                              self.wave_min)
        wave_max = kwargs.get("wave_max",
                              self.wave_max)
        if wave_max < self.wave_min or wave_min > self.wave_max:
            logger.info("   Outside filter range")
            return None
        wave_min = max(wave_min, self.wave_min)
        wave_max = min(wave_max, self.wave_max)
        logger.info("   %.02f .. %.02f um", wave_min, wave_max)

        # bin_width is taken as the minimum dispersion of the trace
        # ..todo: if wcs is given take bin width from cdelt1
        bin_width = kwargs.get("bin_width", None)
        if bin_width is None:
            self._set_dispersion(wave_min, wave_max)
            bin_width = np.abs(self.dlam_per_pix.y).min()
        logger.info("   Bin width %.02g um", bin_width)

        pixscale = from_currsys(self.meta["pixel_scale"], self.cmds)

        # Temporary solution to get slit length
        xi_min = kwargs.get("xi_min", None)
        if xi_min is None:
            try:
                xi_min = hdulist[0].header["HIERARCH INS SLIT XIMIN"]
            except KeyError:
                logger.error("xi_min not found")
                return None
        xi_max = kwargs.get("xi_max", None)
        if xi_max is None:
            try:
                xi_max = hdulist[0].header["HIERARCH INS SLIT XIMAX"]
            except KeyError:
                logger.error("xi_max not found")
                return None

        if wcs is None:
            wcs = WCS(naxis=2)
            wcs.wcs.ctype = ["WAVE", "LINEAR"]
            wcs.wcs.cunit = ["um", "arcsec"]
            wcs.wcs.crpix = [1, 1]
            wcs.wcs.cdelt = [bin_width, pixscale]  # PIXSCALE

        # crval set to wave_min to catch explicitely set value
        wcs.wcs.crval = [wave_min, xi_min]   # XIMIN

        nlam = int((wave_max - wave_min) / bin_width) + 1
        nxi = int((xi_max - xi_min) / pixscale) + 1

        # Create interpolation functions if not provided
        if interps is None:
            logger.info("Computing image interpolations")
            interps = make_image_interpolations(hdulist, kx=1, ky=1)

        # Create Xi, Lam images (do I need Iarr and Jarr or can I build Xi, Lam directly?)
        Iarr, Jarr = np.meshgrid(np.arange(nlam, dtype=np.float32),
                                 np.arange(nxi, dtype=np.float32))
        Lam, Xi = wcs.all_pix2world(Iarr, Jarr, 0)

        # Make sure that we do have microns
        Lam = Lam * u.Unit(wcs.wcs.cunit[0]).to(u.um)

        # Convert Xi, Lam to focal plane units
        Xarr = self.xilam2x(Xi, Lam)
        Yarr = self.xilam2y(Xi, Lam)

        rect_spec = np.zeros_like(Xarr, dtype=np.float32)

        ihdu = 0
        # TODO: make this more iteratory
        for hdu in hdulist:
            if not isinstance(hdu, fits.ImageHDU):
                continue

            wcs_fp = WCS(hdu.header, key="D")
            n_x = hdu.header["NAXIS1"]
            n_y = hdu.header["NAXIS2"]
            iarr, jarr = wcs_fp.all_world2pix(Xarr, Yarr, 0)
            mask = (iarr > 0) * (iarr < n_x) * (jarr > 0) * (jarr < n_y)
            if np.any(mask):
                specpart = interps[ihdu](jarr, iarr, grid=False)
                rect_spec += specpart * mask

            ihdu += 1

        header = wcs.to_header()
        header["EXTNAME"] = self.trace_id
        return fits.ImageHDU(data=rect_spec, header=header)

    def footprint(self, wave_min=None, wave_max=None, xi_min=None, xi_max=None):
        """
        Return corners of rectangle enclosing spectral trace.

        Parameters
        ----------
        wave_min, wave_max : float [um], Quantity
            Minimum and maximum wavelength to compute the footprint on.
            If `None`, use the full range that spectral trace is defined on.
            Float values are interpreted as microns.
        xi_min, xi_max : float [arcsec], Quantity
            Minimum and maximum slit position on the sky.
            If `None`, use the full range that spectral trace is defined on.
            Float values are interpreted as arcsec.
        """
        # Define the wavelength range of the footprint. This is a compromise
        # between the requested range (by method args) and the definition
        # range of the spectral trace
        # This is only relevant if the trace is given by a table of reference
        # points. Otherwise (METIS LMS!) we assume that the range is valid.
        if ("wave_colname" in self.meta and
                self.meta["wave_colname"] in self.table.colnames):
            # Here, the parameters are obtained from a table of reference points
            wave_unit = self.table[self.meta["wave_colname"]].unit
            wave_val = quantify(self.table[self.meta["wave_colname"]].data,
                                wave_unit)

            if wave_min is None:
                wave_min = np.min(wave_val)
            if wave_max is None:
                wave_max = np.max(wave_val)

            wave_min = quantify(wave_min, u.um)
            wave_max = quantify(wave_max, u.um)

            # Requested wavelenth range is entirely outside definition range:
            # no footprint
            if wave_min > np.max(wave_val) or wave_max < np.min(wave_val):
                return None, None

            # Restrict to overlap of requested range and definition range
            wave_min = max(wave_min, np.min(wave_val)).value
            wave_max = min(wave_max, np.max(wave_val)).value

            # Define the slit range of the footprint. This is a compromise
            # between the requested range (by method args) and the definition
            # range of the spectral trace
            try:
                xi_unit = self.table[self.meta["s_colname"]].unit
            except KeyError:
                xi_unit = u.arcsec

            xi_val = quantify(self.table[self.meta["s_colname"]].data,
                              xi_unit)

            if xi_min is None:
                xi_min = np.min(xi_val)
            if xi_max is None:
                xi_max = np.max(xi_val)

            xi_min = quantify(xi_min, u.arcsec)
            xi_max = quantify(xi_max, u.arcsec)

            # Requested slit range is entirely outside definition range:
            # no footprint
            if xi_min > np.max(xi_val) or xi_max < np.min(xi_val):
                return None, None

            # Restrict to overlap of requested range and definition range
            xi_min = max(xi_min, np.min(xi_val)).value
            xi_max = min(xi_max, np.max(xi_val)).value

        # Map the edges of xi/lam to the focal plance
        n_edge = 512
        wave_edge = np.concatenate((np.linspace(wave_min, wave_max, n_edge),
                                    [wave_max] * n_edge,
                                    np.linspace(wave_min, wave_max, n_edge),
                                    [wave_min] * n_edge))
        xi_edge = np.concatenate(([xi_min] * n_edge,
                                  np.linspace(xi_min, xi_max, n_edge),
                                  [xi_max] * n_edge,
                                  np.linspace(xi_min, xi_max, n_edge)))

        x_edge = self.xilam2x(xi_edge, wave_edge)
        y_edge = self.xilam2y(xi_edge, wave_edge)

        return ([x_edge.min(), x_edge.max(), x_edge.max(), x_edge.min()],
                [y_edge.min(), y_edge.min(), y_edge.max(), y_edge.max()])

    def plot(self, wave_min=None, wave_max=None, xi_min=None, xi_max=None, *,
             c="r", axes=None, plot_footprint=True, plot_wave=True,
             plot_ctrlpnts=True, plot_outline=False, plot_trace_id=False):
        """Plot control points (and/or footprint) of the SpectralTrace.

        Parameters
        ----------
        wave_min : float, optional
            Minimum wavelength, if any.
        wave_max : float, optional
            Maximum wavelength, if any.
        xi_min : float, optional
            Minimum slit, if any.
        xi_max : float, optional
            Maximum slit, if any.
        c : str, optional
            Colour, any valid matplotlib colour string. The default is "r".
        axes : matplotlib axes, optional
            The axes object to use for the plot. If None (default), a new
            figure with one axes will be created.

        Returns
        -------
        axes : matplotlib axes
            The axes object containing the plot.

        Other Parameters
        ----------------
        plot_footprint : bool, optional
            Plot a rectangle encompassing all control points, which may be
            larger than the area actually covered by the trace, if the trace is
            not exactly perpendicular to the detector. The default is True.
        plot_wave : bool, optional
            Annotate the wavelength points. The default is True.
        plot_ctrlpnts : bool, optional
            Plot the individual control points as makers. The default is True.
        plot_outline : bool, optional
            Plot the smallest tetragon encompassing all control points.
            The default is False.
        plot_trace_id : bool, optional
            Write the trace ID in the middle of the trace.
            The default is False.
        """
        if axes is None:
            _, axes = figure_factory()

        # Footprint (rectangle enclosing the trace)
        xlim, ylim = self.footprint(wave_min=wave_min, wave_max=wave_max)
        if xlim is None:
            return axes

        if plot_footprint:
            axes.plot(list(close_loop(xlim)), list(close_loop(ylim)))

        # for convenience...
        xname = self.meta["x_colname"]
        yname = self.meta["y_colname"]
        wname = self.meta["wave_colname"]
        sname = self.meta["s_colname"]

        # Control points
        waves = self.table[wname]
        if wave_min is None:
            wave_min = waves.min()
        if wave_max is None:
            wave_max = waves.max()
        xis = self.table[sname]
        if xi_min is None:
            xi_min = xis.min()
        if xi_max is None:
            xi_max = xis.max()

        mask = ((waves >= wave_min) & (waves <= wave_max)
                & (xis >= xi_min) & (xis <= xi_max))
        if sum(mask) <= 2:
            return axes

        w = waves[mask]

        x = self.table[xname][mask]
        y = self.table[yname][mask]
        if plot_ctrlpnts:
            axes.plot(x, y, "o", c=c)

        if plot_outline:
            blue_end = self.table[mask][w == w.min()]
            red_end = self.table[mask][w == w.max()]
            blue_end.sort(sname)
            red_end.sort(sname)
            corners = vstack([blue_end[[0, -1]][xname, yname],
                              red_end[[-1, 0]][xname, yname],
                              blue_end[0][xname, yname]])
            axes.plot(corners[xname], corners[yname], c=c)

        if plot_trace_id:
            axes.text(corners[xname][:-1].mean(), corners[yname][:-1].mean(),
                      self.trace_id, c=c, rotation="vertical",
                      ha="center", va="center")

        for wave in np.unique(w):
            xx = x[w == wave]
            xx.sort()
            dx = xx[-1] - xx[-2]

            if plot_wave:
                axes.text(x[w == wave].max() + 0.5 * dx,
                          y[w == wave].mean(),
                          str(wave), va="center", ha="left")

        axes.set_aspect("equal")
        return axes

    def _set_dispersion(self, wave_min, wave_max, pixsize=None):
        """Compute of dispersion dlam_per_pix along xi=0."""
        # ..todo: This may have to be generalised - xi=0 is at the centre
        # of METIS slits and the short MICADO slit.

        xi = np.array([0] * 1001)
        lam = np.linspace(wave_min, wave_max, 1001)
        x_mm = self.xilam2x(xi, lam)
        y_mm = self.xilam2y(xi, lam)
        if self.dispersion_axis == "x":
            dlam_grad = self.xy2lam.gradient()[0]  # dlam_by_dx
        else:
            dlam_grad = self.xy2lam.gradient()[1]  # dlam_by_dy
        pixsize = (from_currsys(self.meta["pixel_scale"], self.cmds) /
                   from_currsys(self.meta["plate_scale"], self.cmds))
        self.dlam_per_pix = interp1d(lam,
                                     dlam_grad(x_mm, y_mm) * pixsize,
                                     fill_value="extrapolate")

    def __repr__(self):
        return f"{self.__class__.__name__}({self.table!r}, **{self.meta!r})"

    def __str__(self):
        msg = (f"<SpectralTrace> \"{self.trace_id}\" : "
               f"[{self.wave_min:.4f}, {self.wave_max:.4f}]um : "
               f"Ext {self.meta['extension_id']} : "
               f"Aperture {self.meta['aperture_id']} : "
               f"ImagePlane {self.meta['image_plane_id']}")
        return msg


class XiLamImage():
    """
    Class to compute a rectified 2D spectrum.

    The class produces and holds an image of xi (relative position along
    the spatial slit direction) and wavelength lambda.

    Parameters
    ----------
    fov : FieldOfView
    dlam_per_pix : a 1-D interpolation function from wavelength (in um) to
          dispersion (in um/pixel); alternatively a number giving an average
          dispersion
    """

    def __init__(self, fov, dlam_per_pix):
        # ..todo: we assume that we always have a cube. We use SpecCADO's
        #         add_cube_layer method
        cube_wcs = WCS(fov.cube.header, key=" ")
        wcs_lam = cube_wcs.sub([3])

        d_xi = fov.cube.header["CDELT1"]
        d_xi *= u.Unit(fov.cube.header["CUNIT1"]).to(u.arcsec)
        d_eta = fov.cube.header["CDELT2"]
        d_eta *= u.Unit(fov.cube.header["CUNIT2"]).to(u.arcsec)
        d_lam = fov.cube.header["CDELT3"]
        d_lam *= u.Unit(fov.cube.header["CUNIT3"]).to(u.um)

        # This is based on the cube shape and assumes that the cube's spatial
        # dimensions are set by the slit aperture
        (n_lam, n_eta, n_xi) = fov.cube.data.shape

        # arrays of cube coordinates
        cube_xi = d_xi * np.arange(n_xi) + fov.meta["xi_min"].value
        cube_eta = d_eta * (np.arange(n_eta) - (n_eta - 1) / 2)
        cube_lam = wcs_lam.all_pix2world(np.arange(n_lam), 1)[0]
        cube_lam *= u.Unit(wcs_lam.wcs.cunit[0]).to(u.um)

        # Initialise the array to hold the xi-lambda image
        self.image = np.zeros((n_xi, n_lam), dtype=np.float32)
        self.lam = cube_lam
        try:
            dlam_per_pix_val = dlam_per_pix(np.asarray(self.lam))
        except TypeError:
            dlam_per_pix_val = dlam_per_pix
            logger.warning("Using scalar dlam_per_pix = %.2g",
                           dlam_per_pix_val)

        for i, eta in enumerate(cube_eta):
            lam0 = self.lam + dlam_per_pix_val * eta / d_eta

            # lam0 is the target wavelength. We need to check that this
            # overlaps with the wavelength range covered by the cube
            if lam0.min() < cube_lam.max() and lam0.max() > cube_lam.min():
                plane = fov.cube.data[:, i, :].T
                plane_interp = RectBivariateSpline(cube_xi, cube_lam, plane,
                                                   kx=1, ky=1)
                self.image += plane_interp(cube_xi, lam0)

        self.image *= d_eta     # ph/s/um/arcsec2 --> ph/s/um/arcsec

        # WCS for the xi-lambda image, i.e. the rectified 2D spectrum
        # Default WCS with xi in arcsec
        self.wcs = WCS(naxis=2)
        self.wcs.wcs.crpix = [1, 1]
        self.wcs.wcs.crval = [self.lam[0], fov.meta["xi_min"].value]
        self.wcs.wcs.pc = [[1, 0], [0, 1]]
        self.wcs.wcs.cdelt = [d_lam, d_xi]
        self.wcs.wcs.ctype = ["LINEAR", "LINEAR"]
        self.wcs.wcs.cname = ["WAVELEN", "SLITPOS"]
        self.wcs.wcs.cunit = ["um", "arcsec"]

        # Alternative: xi = [0, 1], dimensionless
        self.wcsa = WCS(naxis=2)
        self.wcsa.wcs.crpix = [1, 1]
        self.wcsa.wcs.crval = [self.lam[0], 0]
        self.wcsa.wcs.pc = [[1, 0], [0, 1]]
        self.wcsa.wcs.cdelt = [d_lam, 1./n_xi]
        self.wcsa.wcs.ctype = ["LINEAR", "LINEAR"]
        self.wcsa.wcs.cname = ["WAVELEN", "SLITPOS"]
        self.wcs.wcs.cunit = ["um", ""]

        self.xi = self.wcs.all_pix2world(self.lam[0], np.arange(n_xi), 0)[1]
        self.npix_xi = n_xi
        self.npix_lam = n_lam
        # ..todo: cubic spline introduces negative values, linear does not.
        #  Alternative might be to cubic-spline interpolate on sqrt(image),
        #  with subsequent squaring of the result. This would require
        #  wrapping RectBivariateSpline in a new (sub)class.
        spline_order = (1, 1)
        self.interp = RectBivariateSpline(self.xi, self.lam, self.image,
                                          kx=spline_order[0],
                                          ky=spline_order[1])


class Transform2D():
    """
    2-dimensional polynomial transform.

    The class is instantiated from a m x n matrix A that contains the
    coefficients of the polynomial. Along rows, the power of x increases;
    along columns, the power of y increases, such that A[j, i] is
    the coefficient of x^i y^j.

    The functions `pretransform_x` and `pretransform_y` can be used to
    transform the input variables before the matrix is applied. The function
    `posttransform` can be applied to the output after application of the
    matrix.

    In Scopesim, a usecase for the pre- and post-transform functions is the
    METIS LMS, where the matrices are applied to phases while Scopesim
    operates on wavelengths. The functions to pass are `lam2phase` and
    `phase2lam`.

    Parameters
    ----------
    matrix : np.array
        matrix of polynomial coefficients
    pretransform_x : function, tuple
    pretransform_y : function, tuple
        If not None, the function is applied to the input
        variable `x` or `y` before the actual 2D transform is computed
    posttransform : function, tuple
        If not None, the function is applied to the output variable
        after the 2D transform is computed

    When passed as a tuple, the first element is the function itself,
    the second element is a dictionary of arguments to the function.
    Example:
    ```
    def rescale(x, scale=1.):
        return x * scale
    pretransform_x = (rescale, {"scale": 0.5})
    ```
    """

    def __init__(self, matrix, pretransform_x=None,
                 pretransform_y=None, posttransform=None):
        self.matrix = np.asarray(matrix)
        self.ny, self.nx = self.matrix.shape
        self.pretransform_x = self._repackage(pretransform_x)
        self.pretransform_y = self._repackage(pretransform_y)
        self.posttransform = self._repackage(posttransform)

    def _repackage(self, trafo):
        """Make sure `trafo` is a tuple."""
        if trafo is not None and not isinstance(trafo, tuple):
            trafo = (trafo, {})
        return trafo

    def __call__(self, x, y, grid=False, **kwargs):
        """
        Apply the polynomial transform.

        The transformation is a polynomial based on the simple
        monomials x^i y^j. When grid=True, the transform is applied to the grid
        formed by the tensor product of the vectors x and y. When grid=False,
        the vectors of x and y define the components of a number of points (in
        this case, x and y must be of the same length).

        Functions `pretransform_x`, `pretransform_y` and `posttransform`
        can be supplied to override the instance values.

        Parameters
        ----------
        x, y : np.array
            x and y values to transform
        grid : boolean
            If true, return result for all pairs of components in x and y.

        Return
        ------
        When grid=True, a matrix with results for all pairings of components
        in x and y. When grid=False, a vector. In this case, x and y must
        have the same length.
        """
        if "pretransform_x" in kwargs:
            self.pretransform_x = self._repackage(kwargs["pretransform_x"])
        if "pretransform_y" in kwargs:
            self.pretransform_y = self._repackage(kwargs["pretransform_y"])
        if "posttransform" in kwargs:
            self.posttransform = self._repackage(kwargs["posttransform"])

        x = np.array(x)
        y = np.array(y)
        orig_shape = x.shape

        if not grid and x.shape != y.shape:
            raise ValueError("x and y must have the same length when grid "
                             "is False")

        # Apply pre transforms
        if self.pretransform_x is not None:
            x = self.pretransform_x[0](x, **self.pretransform_x[1])
        if self.pretransform_y is not None:
            y = self.pretransform_y[0](y, **self.pretransform_y[1])

        xvec = power_vector(x.flatten(), self.nx - 1)
        yvec = power_vector(y.flatten(), self.ny - 1)

        temp = self.matrix @ xvec

        if grid:
            result = yvec.T @ temp
        else:
            # Compute the scalar product of each column in yvec with the
            # corresponding column in temp. This gives the diagonal of the
            # expression in the "grid" branch.
            result = (yvec * temp).sum(axis=0)
            if not orig_shape:
                result = np.float32(result)
            else:
                result = result.reshape(orig_shape)

        # Apply posttransform
        if self.posttransform is not None:
            result = self.posttransform[0](result, **self.posttransform[1])

        return result

    @classmethod
    def fit(cls, xin, yin, xout, degree=4):
        """Determine polynomial fits."""
        pinit = Polynomial2D(degree=degree)
        fitter = fitting.LinearLSQFitter()
        fit = fitter(pinit, xin, yin, xout)
        return Transform2D(fit2matrix(fit))

    def gradient(self):
        """Compute the gradient of a 2d polynomial transformation."""
        mat = self.matrix

        dmat_x = (mat * np.arange(self.nx))[:, 1:]
        dmat_y = (mat.T * np.arange(self.ny)).T[1:, :]

        return Transform2D(dmat_x), Transform2D(dmat_y)


def fit2matrix(fit):
    """
    Return coefficients from a polynomial fit as a matrix.

    The Polynomial2D fits of degree n have coefficients for all i, j
    with i + j <= n.
    How would one rearrange those?
    """
    coeffs = dict(zip(fit.param_names, fit.parameters))
    deg = fit.degree
    mat = np.zeros((deg + 1, deg + 1), dtype=np.float32)
    for i in range(deg + 1):
        for j in range(deg + 1):
            try:
                mat[j, i] = coeffs[f"c{i}_{j}"]
            except KeyError:
                pass
    return mat


# ..todo: should the next three functions be combined and return a dictionary of fits?
def xilam2xy_fit(layout, params):
    """
    Determine polynomial fits of FPA position.

    Fits are of degree 4 as a function of slit position and wavelength.
    """
    xi_arr = layout[params["s_colname"]]
    lam_arr = layout[params["wave_colname"]]
    x_arr = layout[params["x_colname"]]
    y_arr = layout[params["y_colname"]]

    # Filter the lists: remove any points with x==0
    # ..todo: this may not be necessary after sanitising the table
    # good = x != 0
    # xi = xi[good]
    # lam = lam[good]
    # x = x[good]
    # y = y[good]

    # compute the fits
    pinit_x = Polynomial2D(degree=4)
    pinit_y = Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xilam2x = fitter(pinit_x, xi_arr, lam_arr, x_arr)
    xilam2y = fitter(pinit_y, xi_arr, lam_arr, y_arr)

    return xilam2x, xilam2y


def xy2xilam_fit(layout, params):
    """
    Determine polynomial fits of wavelength/slit position.

    Fits are of degree 4 as a function of focal plane position
    """
    xi_arr = layout[params["s_colname"]]
    lam_arr = layout[params["wave_colname"]]
    x_arr = layout[params["x_colname"]]
    y_arr = layout[params["y_colname"]]

    pinit_xi = Polynomial2D(degree=4)
    pinit_lam = Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xy2xi = fitter(pinit_xi, x_arr, y_arr, xi_arr)
    xy2lam = fitter(pinit_lam, x_arr, y_arr, lam_arr)

    return xy2xi, xy2lam


def _xiy2xlam_fit(layout, params):
    """Determine polynomial fits of wavelength/slit position.

    Fits are of degree 4 as a function of focal plane position
    """
    # These are helper functions to allow fitting of left/right edges
    # for the purpose of checking whether a trace is on a chip or not.

    xi_arr = layout[params["s_colname"]]
    lam_arr = layout[params["wave_colname"]]
    x_arr = layout[params["x_colname"]]
    y_arr = layout[params["y_colname"]]

    pinit_x = Polynomial2D(degree=4)
    pinit_lam = Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xiy2x = fitter(pinit_x, xi_arr, y_arr, x_arr)
    xiy2lam = fitter(pinit_lam, xi_arr, y_arr, lam_arr)
    return xiy2x, xiy2lam


def make_image_interpolations(hdulist, **kwargs):
    """Create 2D interpolation functions for images."""
    interps = []
    for hdu in hdulist:
        if isinstance(hdu, fits.ImageHDU):
            interps.append(
                RectBivariateSpline(np.arange(hdu.header["NAXIS1"]),
                                    np.arange(hdu.header["NAXIS2"]),
                                    hdu.data, **kwargs)
            )
    return interps

# ..todo: Check whether the following functions are actually used


def rolling_median(x, n):
    """Calculate the rolling median of a sequence for +/- n entries."""
    y = [np.median(x[max(0, i-n):min(len(x), i+n+1)]) for i in range(len(x))]
    return np.array(y)


def fill_zeros(x):
    """Fill in zeros in a sequence with the previous non-zero number."""
    for i in range(1, len(x)):
        if x[i] == 0:
            x[i] = x[i-1]
    return x


def get_affine_parameters(coords):
    """
    Return rotation and shear for each MTC point along a ``SpectralTrace``.

    .. note: Restrictions of this method:

       * only uses the left most coordinates for the shear
       * rotation angle is calculated using the trace extremes

    Parameters
    ----------
    coords : dict of 2D arrays
        Each dict entry ["x", "y", "s"] contains a [N, M] 2D array of
        coordinates, where:

        * N is the number of points along the slit (e.g. ~5), and
        * M is the number of positions along the trace (e.g. >100)

    Returns
    -------
    rotations : array
        [deg] Rotation angles for M positions along the Trace
    shears : array
        [deg] Shear angles for M positions along the Trace

    """
    rad2deg = 180 / np.pi
    dxs = coords["x"][-1, :] - coords["x"][0, :]
    dys = coords["y"][-1, :] - coords["y"][0, :]
    rotations = np.arctan2(dys, dxs) * rad2deg

    dxs = np.diff(coords["x"], axis=1)
    dys = np.diff(coords["y"], axis=1)
    shears = np.array([np.arctan2(dys[i], dxs[i])
                      for i in range(dxs.shape[0])])
    shears = np.array(list(shears.T) + [shears.T[-1]]).T
    shears = (np.average(shears, axis=0) * rad2deg) - (90 + rotations)

    return rotations, shears
