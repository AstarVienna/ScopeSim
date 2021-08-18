"""
This module contains the definition of the `SpectralTrace` class.
"""
from copy import deepcopy

import numpy as np

from scipy.interpolate import RectBivariateSpline
from matplotlib import pyplot as plt

from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

from ..optics import spectral_trace_utils as spt_utils
from ..optics import image_plane_utils as imp_utils
from ..optics.monochromatic_trace_curve import MonochromeTraceCurve
from ..utils import interp2, check_keys, from_currsys, quantify


class SpectralTrace:
    '''Definition of one spectral trace

    A SpectralTrace describes the mapping of spectral slit coordinates
    to the focal plane. The class reads an order layout and fits several
    functions to describe the geometry of the trace.

    Slit coordinates are:
    - xi : spatial position along the slit [arcsec]
    - lam : Wavelength [um]
    Focal plane coordinates are:
    - x, y : [mm]
    '''

    def __init__(self, trace_tbl, **kwargs):
        self.meta = {"x_colname": "x",
                     "y_colname": "y",
                     "s_colname": "s",
                     "wave_colname": "wavelength",
                     "col_number_start": 0,
                     "dwave": 0.002,
                     "aperture_id": 0,
                     "image_plane_id": 0,
                     "extension_id": 2,
                     "invalid_value": None,
                     "spline_order": 4,
                     "pixel_size": None,
                     "description": "<no description>"}
        self.meta.update(kwargs)

        if isinstance(trace_tbl, (fits.BinTableHDU, fits.TableHDU)):
            self.table = Table.read(trace_tbl)
            try:
                self.meta["description"] = trace_tbl.header['EXTNAME']
            except KeyError:
                pass
        elif isinstance(trace_tbl, Table):
            self.table = trace_tbl
        else:
            raise ValueError("trace_tbl must be one of (fits.BinTableHDU, "
                             "fits.TableHDU, astropy.Table): {}"
                             "".format(type(trace_tbl)))

        if self.meta["invalid_value"] is not None:
            self.table = spt_utils.sanitize_table(
                self.table,
                invalid_value=self.meta["invalid_value"],
                wave_colname=self.meta["wave_colname"],
                x_colname=self.meta["x_colname"],
                y_colname=self.meta["y_colname"],
                spline_order=self.meta["spline_order"],
                ext_id=self.meta["extension_id"])
        # ..todo: Should that be made np.unique?
        self.waves = self.table[self.meta["wave_colname"]]
        #..todo: this turns out as 1 in simplified layout. Remove?
        self.n_traces = len([col for col in self.table.colnames
                             if self.meta["y_colname"] in col])

        self.wave_min = quantify(np.min(self.waves), u.um).value
        self.wave_max = quantify(np.max(self.waves), u.um).value

        # ..todo: This bunch of variables can probably go.
        self._disp = None           # [mm/um] spectral dispersion distance
        self._waverange = None
        self._wave_bin_edges = None
        self._wave_bin_centers = None
        self._curves = None

        # Interpolation functions
        self.xy2xi, self.xy2lam = spt_utils.xy2xilam_fit(self.table,
                                                         self.meta)
        self.xilam2x, self.xilam2y = spt_utils.xilam2xy_fit(self.table,
                                                            self.meta)
        self._xiy2x, self._xiy2lam = spt_utils._xiy2xlam_fit(self.table,
                                                             self.meta)

    def map_spectra_to_focal_plane(self, fov):
        """
        Apply the spectral trace mapping to a spectral cube

        The cube is contained in a FieldOfView object, which also has
        world coordinate systems for the Source (sky coordinates and
        wavelengths) and for the focal plane.
        The method returns a section of the fov image along with info on
        where this image lies in the focal plane.
        """
        # Initialise the image based on the footprint of the spectral
        # trace and the focal plane WCS
        wave_min = fov.meta['wave_min']
        wave_max = fov.meta['wave_max']
        xi_min = fov.meta['xi_min']
        xi_max = fov.meta['xi_max']
        xlim, ylim = self.footprint(wave_min=wave_min, wave_max=wave_max,
                                    xi_min=xi_min, xi_max=xi_max)

        if xlim is None:
            return None

        # WCSD from the FieldOfView
        wcsd = WCS(fov.header, key='D')
        naxis1, naxis2 = fov.header['NAXIS1'], fov.header['NAXIS2']
        xpix, ypix = wcsd.all_world2pix(xlim, ylim, 0) ## oder 1?
        xmin = np.floor(xpix.min()).astype(int)
        xmax = np.ceil(xpix.max()).astype(int)
        ymin = np.floor(ypix.min()).astype(int)
        ymax = np.ceil(ypix.max()).astype(int)

        # Check if spectral trace footprint is outside FoV
        if xmax < 0 or xmin > naxis1 or ymax < 0 or ymin > naxis2:
            return None, xmin, xmax, ymin, ymax

        # Only work on parts within the FoV
        xmin = max(xmin, 0)
        xmax = min(xmax, naxis1)
        ymin = max(ymin, 0)
        ymax = min(ymax, naxis2)

        # Update header for the subimage
        fov.header['CRPIX1'] -= xmin
        fov.header['CRPIX2'] -= ymin
        fov.header['CRPIX1D'] -= xmin
        fov.header['CRPIX2D'] -= ymin
        fov.header['NAXIS1'] = xmax - xmin
        fov.header['NAXIS2'] = ymax - ymin

        # initialise the subimage
        image = np.zeros((fov.header['NAXIS2'], fov.header['NAXIS1']), dtype=np.float32)

        # WCSD from the FieldOfView - this is now for the subimage
        wcsd = WCS(fov.header, key='D')

        xmin_um, ymin_um = wcsd.all_pix2world(xmin, ymin, 0)
        # xmax_um, ymax_um = wcsd.all_pix2world(xmax, ymax, 0)
        pixsize = fov.header['CDELT1D']
        pix_edges_x = xmin_um + np.arange(0, fov.header['NAXIS1'] + 1) * pixsize
        pix_edges_y = ymin_um + np.arange(0, fov.header['NAXIS2'] + 1) * pixsize

        # Step 512 lines at a time ..todo: necessary? Do full image first
        # ..todo: this appears to be the old version with dispersion in the y direction
        im_j1 = 0
        im_j2 = fov.header['NAXIS2']

        ymin = pix_edges_y[im_j1]
        ymax = pix_edges_y[im_j2]

        # wavelength range for the image slice wave_min, wave_max
        # ..todo: in speccado loop over 512 rows, wavelength limits for each slice
        xcen = self.xilam2x((xi_min + xi_max) / 2, (wave_min + wave_max) / 2)

        # wavelength step per detector pixel at centre of slice
        # ..todo: I try the average here, but is that correct? The pixel
        # edges may not correspond precisely to wave limits
        dlam_per_pix = (wave_max - wave_min) / fov.header['NAXIS2']
        try:
            xilam = XiLamImage(fov, dlam_per_pix)
        except ValueError:
            print(" ---> " + self.description + "gave ValueError")

        npix_xi, npix_lam = xilam.npix_xi, xilam.npix_lam

        lam = xilam.lam
        xi = xilam.xi

        image_interp = xilam.interp
        ## ..todo: continue simulation.py l.208


        fov._image = image

    def get_max_dispersion(self, **kwargs):
        '''Get the maximum dispersion in a spectral trace

        This is a wrapper for the function in `spectral_trace_utils`.
        '''
        params = {}
        params.update(self.meta)
        params["wave_min"] = self.wave_min
        params["wave_max"] = self.wave_max
        params.update(kwargs)

        disp, waverange = spt_utils.get_max_dispersion(
            self, **params)  # dwave is passed from kwargs
        self._disp = disp
        self._waverange = waverange

        return self._disp, self._waverange

    def get_pixel_wavelength_edges(self, pixel_size):
        """Returns the wavelengths at the edge of pixels along a trace"""
        self.meta["pixel_size"] = pixel_size
        if self._disp is None or self._waverange is None:
            self.get_max_dispersion()

        um_per_pix = pixel_size / self._disp  # wavelength range per pixel
        wbe = spt_utils.pixel_wavelength_edges(um_per_pix, self._waverange,
                                               self.wave_min, self.wave_max)
        self._wave_bin_edges = wbe
        self._wave_bin_centers = 0.5 * (wbe[:-1] + wbe[1:])

        return self._wave_bin_edges

    @property
    def wave_edges(self):
        pixel_size = self.meta["pixel_size"]
        if self._wave_bin_edges is None and pixel_size is not None:
            self.get_pixel_wavelength_edges(pixel_size)
        return self._wave_bin_edges

    @property
    def wave_centers(self):
        pixel_size = self.meta["pixel_size"]
        if self._wave_bin_centers is None and pixel_size is not None:
            self.get_pixel_wavelength_edges(pixel_size)
        return self._wave_bin_centers

    def footprint(self, wave_min=None, wave_max=None, xi_min=None, xi_max=None):
        '''
        Return corners of rectangle enclosing spectral trace

        Parameters
        ----------
        wave_min, wave_max : float [um], Quantity
            Minimum and maximum wavelength to compute the footprint on.
            If `None`, use the full range that the spectral trace is defined on.
            Float values are interpreted as microns.
        xi_min, xi_max : float [arcsec], Quantity
            Minimum and maximum slit position on the sky.
            If `None`, use the full range that the spectral trace is defined on.
            Float values are interpreted as arcsec.
        '''
        ## Define the wavelength range of the footprint. This is a compromise
        ## between the requested range (by method args) and the definition
        ## range of the spectral trace
        try:
            wave_unit = self.table[self.meta['wave_colname']].unit
        except KeyError:
            wave_unit = u.um

        wave_val = quantify(self.table[self.meta['wave_colname']].data,
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

        ## Define the slit range of the footprint. This is a compromise
        ## between the requested range (by method args) and the definition
        ## range of the spectral trace
        try:
            xi_unit = self.table[self.meta['s_colname']].unit
        except KeyError:
            xi_unit = u.arcsec

        xi_val = quantify(self.table[self.meta['s_colname']].data,
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

        return ([np.min(x_edge), np.max(x_edge), np.max(x_edge), np.min(x_edge)],
                [np.min(y_edge), np.min(y_edge), np.max(y_edge), np.max(y_edge)])

    def get_trace_curves(self, pixel_size, wave_min=None, wave_max=None,
                         xy_edges=None):
        """Returns a list of MonochromeTraceCurves for projecting apertures

        Parameters
        ----------
        pixel_size : float
            [mm]
        wave_min, wave_max : float
            [um]
        xy_edges : list of floats
            [mm] Borders of the usable region. Ignore trace curves outside
        """
        wave_min = self.wave_min if wave_min is None else wave_min
        wave_max = self.wave_max if wave_max is None else wave_max
        if self._wave_bin_edges is None:
            self.get_pixel_wavelength_edges(pixel_size)

        mask = (self._wave_bin_edges >= wave_min) * \
               (self._wave_bin_edges <= wave_max)

        if sum(mask) == 0:
            self._curves = []
            return self._curves

        wave_edges = self._wave_bin_edges[mask]
        wave_cens = 0.5 * (wave_edges[:-1] + wave_edges[1:])

        # If we want to add polynomial descriptions of the trace curves,
        # here is where it will go
        #
        # We need another construct other than self.table plus another function
        # which creates the x,y,s coordinates for each MonochromeTraceCurve

        n = self.n_traces
        k = self.meta["col_number_start"]
        coords = {z: None for z in "xys"}
        for z in "xys":
            cols = [self.meta[z+"_colname"] + str(ii) for ii in range(k, n+k)]
            coords[z] = np.array([interp2(wave_cens, self.waves,
                                          self.table[col]) for col in cols])

        if xy_edges is not None:
            # ..todo:: not perfect. Ignores traces where two points are ouside,
            #          but the line between them crosses into the xy_region.
            mask = np.any((coords["x"] >= xy_edges["x_min"]) *
                          (coords["x"] <= xy_edges["x_max"]) *
                          (coords["y"] >= xy_edges["y_min"]) *
                          (coords["y"] <= xy_edges["y_max"]), axis=0)
        else:
            mask = [True] * len(wave_cens)

        # ..todo: Don't like this - fix it!
        if sum(mask) == 0:
            self._curves = []
            return self._curves

        rotation, shear = spt_utils.get_affine_parameters(coords)
        self._curves = [MonochromeTraceCurve(x=coords["x"][:, ii],
                                             y=coords["y"][:, ii],
                                             s=coords["s"][:, ii],
                                             wave_min=wave_edges[ii],
                                             wave_max=wave_edges[ii+1],
                                             pixel_size=pixel_size,
                                             rotation=rotation[ii],
                                             shear=shear[ii])
                        for ii in range(len(wave_cens)) if mask[ii]]

        return self._curves

    def get_curve_headers(self, pixel_size, wave_min=None, wave_max=None,
                          detector_edges=None):
        """Collect all the headers from the list of MonochromaticTraceCurves"""
        if self._curves is None:
            self.get_trace_curves(pixel_size, wave_min, wave_max,
                                  xy_edges=detector_edges)
        headers = [curve.get_header(pixel_size) for curve in self._curves]
        return headers

    def fov_headers(self, sky_header, **kwargs):
        check_keys(kwargs, ["wave_min", "wave_max",
                            "pixel_scale", "plate_scale"], "error")
        kwargs = from_currsys(kwargs)
        wave_min = quantify(kwargs["wave_min"], u.um).value
        wave_max = quantify(kwargs["wave_max"], u.um).value
        pixel_D_scale = sky_header["CDELT1"] / (kwargs["plate_scale"] / 3600)

        detector_edges = None
        if "det_header" in kwargs and kwargs["det_header"] is not None:
            xdet, ydet = imp_utils.calc_footprint(kwargs["det_header"], "D")
            detector_edges = {"x_min": np.min(xdet), "x_max": np.max(xdet),
                              "y_min": np.min(ydet), "y_max": np.max(ydet)}

        if sky_header["APERTURE"] != self.meta["aperture_id"]:
            fov_hdrs = []
        elif wave_min > self.wave_max or wave_max < self.wave_min:
            fov_hdrs = []
        else:
            pixel_size = kwargs["pixel_scale"] / kwargs["plate_scale"]
            curve_hdrs = self.get_curve_headers(pixel_size, wave_min, wave_max,
                                                detector_edges=detector_edges)
            if len(curve_hdrs) > 0:
                print("Generated {} headers from {}".format(len(curve_hdrs),
                                                            self.__repr__()))

            for mtc_hdr in curve_hdrs:
                mtc_hdr["EXT"] = self.meta["extension_id"]
                mtc_hdr["APERTURE"] = self.meta["aperture_id"]
                mtc_hdr["IMGPLANE"] = self.meta["image_plane_id"]
                mtc_hdr["CDELT1D"] = pixel_D_scale
                mtc_hdr["CDELT2D"] = pixel_D_scale
                # mtc_hdr["CRPIX2D"] = sky_header["NAXIS2"] * 0.5
                # ..todo:: assumption here is that they are on the same pixel scale - bad assumption!
                mtc_hdr["NAXIS2"] = sky_header["NAXIS2"]
                mtc_hdr.update(sky_header)

            fov_hdrs = curve_hdrs

        return fov_hdrs

    def plot(self, wave_min=None, wave_max=None, c="r"):
        '''Plot control points of the SpectralTrace'''

        # Footprint (rectangle enclosing the trace)
        xlim, ylim  = self.footprint(wave_min=wave_min, wave_max=wave_max)
        xlim.append(xlim[0])
        ylim.append(ylim[0])
        plt.plot(xlim, ylim)

        # Control points
        waves = self.table[self.meta["wave_colname"]]
        if wave_min is None:
            wave_min = waves.min()
        if wave_max is None:
            wave_max = waves.max()

        mask = (waves >= wave_min) * (waves <= wave_max)
        if sum(mask) > 2:
            w = waves[mask]

            x = self.table[self.meta["x_colname"]][mask]
            y = self.table[self.meta["y_colname"]][mask]
            plt.plot(x, y, 'o', c=c)

            for wave in np.unique(waves):
                xx = x[waves==wave]
                xx.sort()
                dx = xx[-1] - xx[-2]
                plt.text(x[waves==wave].max() + 0.5 * dx,
                         y[waves==wave].mean(),
                         str(wave), va='center', ha='left')


            plt.gca().set_aspect("equal")

    def __repr__(self):
        msg = '<SpectralTrace> "{}" : [{}, {}]um : Ext {} : Aperture {} : ' \
              'ImagePlane {}' \
              ''.format(self.meta["description"],
                        round(self.wave_min, 4), round(self.wave_max, 4),
                        self.meta["extension_id"], self.meta["aperture_id"],
                        self.meta["image_plane_id"])
        return msg



class XiLamImage():
    """
    Class to compute a rectified 2D spectrum

    The class produces and holds an image of xi (relative position along
    the spatial slit direction) and wavelength lambda.
    """

    def __init__(self, fov, dlam_per_pix):

        # Slit dimensions: oversample with respect to detector pixel scale

        # Compute size of xi-lambda image: npix_xi, npix_lam


        # ..todo: we assume that we always have a cube. We use SpecCADO's
        #         add_cube_layer method

        cube_wcs = WCS(fov.cube.header, key=' ')
        wcs_lam = cube_wcs.sub([3])

        d_xi = fov.cube.header['CDELT1']
        d_xi *= u.Unit(fov.cube.header['CUNIT1']).to(u.arcsec)
        d_eta = fov.cube.header['CDELT2']
        d_eta *= u.Unit(fov.cube.header['CUNIT2']).to(u.arcsec)
        d_lam = fov.cube.header['CDELT3']
        d_lam *= u.Unit(fov.cube.header['CUNIT3']).to(u.um)

        # This is based on the cube shape and assumes that the cube's spatial
        # dimensions are set by the slit aperture
        (n_lam, n_eta, n_xi) = fov.cube.data.shape

        # arrays of cube coordinates
        cube_xi = d_xi * np.arange(n_xi) + fov.meta['xi_min'].value
        cube_eta = d_eta * (np.arange(n_eta) - (n_eta - 1) / 2)
        cube_lam = wcs_lam.all_pix2world(np.arange(n_lam), 0)[0]
        cube_lam *= u.Unit(wcs_lam.wcs.cunit[0]).to(u.um)

        # Initialise the array to hold the xi-lambda image
        self.image = np.zeros((n_xi, n_lam), dtype=np.float32)
        self.lam = cube_lam

        for i, eta in enumerate(cube_eta):
            #if abs(eta) > fov.slit_width / 2:   # ..todo: needed?
            #    continue

            lam0 = self.lam + dlam_per_pix * eta / d_eta
            # lam0 is the target wavelength. We need to check that this
            # overlaps with the wavelength range covered by the cube
            if lam0.min() < cube_lam.max() and lam0.max() > cube_lam.min():
                plane = fov.cube.data[:, i, :].T
                plane_interp = RectBivariateSpline(cube_xi, cube_lam, plane)
                self.image += plane_interp(cube_xi, lam0)

        # WCS for the xi-lambda image, i.e. the rectified 2D spectrum
        # Default WCS with xi in arcsec
        self.wcs = WCS(naxis=2)
        self.wcs.wcs.crpix = [1, 1]
        self.wcs.wcs.crval = [self.lam[0], fov.meta['xi_min'].value]
        self.wcs.wcs.pc = [[1, 0], [0, 1]]
        self.wcs.wcs.cdelt = [d_lam, d_xi]
        self.wcs.wcs.ctype = ['LINEAR', 'LINEAR']
        self.wcs.wcs.cname = ['WAVELEN', 'SLITPOS']
        self.wcs.wcs.cunit = ['um', 'arcsec']

        # Alternative: xi = [0, 1], dimensionless
        self.wcsa = WCS(naxis=2)
        self.wcsa.wcs.crpix = [1, 1]
        self.wcsa.wcs.crval = [self.lam[0], 0]
        self.wcsa.wcs.pc = [[1, 0], [0, 1]]
        self.wcsa.wcs.cdelt = [d_lam, 1./n_xi]
        self.wcsa.wcs.ctype = ['LINEAR', 'LINEAR']
        self.wcsa.wcs.cname = ['WAVELEN', 'SLITPOS']
        self.wcs.wcs.cunit = ['um', '']

        self.xi = self.wcs.all_pix2world(self.lam[0], np.arange(n_xi), 0)[1]
        self.npix_xi = n_xi
        self.npix_lam = n_lam
        self.interp = RectBivariateSpline(self.xi, self.lam, self.image)
