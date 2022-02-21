"""
This module contains
   - the definition of the `SpectralTrace` class.
   - the definition of the `XiLamImage` class
   - utility functions for use with spectral traces
"""

import logging

import numpy as np

from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import pyplot as plt

from astropy.table import Table
from astropy.modeling import fitting
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.modeling.models import Polynomial2D

from ..optics import image_plane_utils as imp_utils
from ..utils import deriv_polynomial2d, power_vector, interp2, check_keys,\
    from_currsys, quantify


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
                self.meta["trace_id"] = trace_tbl.header['EXTNAME']
            except KeyError:
                pass
        elif isinstance(trace_tbl, Table):
            self.table = trace_tbl
        else:
            raise ValueError("trace_tbl must be one of (fits.BinTableHDU, "
                             "fits.TableHDU, astropy.Table): {}"
                             "".format(type(trace_tbl)))

        if self.meta["invalid_value"] is not None:
            self.table = sanitize_table(
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

        # Interpolation functions   ..todo: equivalent for LMS
        self.compute_interpolation_functions()

    def fov_grid(self):
        """
        Provide information on the source space volume required by the effect

        Returns
        -------
        A dictionary with entries `wave_min` and `wave_max`.
        Spatial limits are determined by the `ApertureMask` effect
        and are not returned here.
        """
        trace_id = self.meta['trace_id']
        aperture_id = self.meta['aperture_id']
        lam_arr = self.table[self.meta['wave_colname']]

        wave_max = np.max(lam_arr)
        wave_min = np.min(lam_arr)

        return {'wave_min': wave_min, 'wave_max': wave_max,
                'trace_id': trace_id, 'aperture_id': aperture_id}

    def compute_interpolation_functions(self):
        """
        Compute various interpolation functions between slit and focal plane
        """
        x_arr = self.table[self.meta['x_colname']]
        y_arr = self.table[self.meta['y_colname']]
        xi_arr = self.table[self.meta['s_colname']]
        lam_arr = self.table[self.meta['wave_colname']]

        self.xy2xi = Transform2D.fit(x_arr, y_arr, xi_arr)
        self.xy2lam = Transform2D.fit(x_arr, y_arr, lam_arr)
        self.xilam2x = Transform2D.fit(xi_arr, lam_arr, x_arr)
        self.xilam2y = Transform2D.fit(xi_arr, lam_arr, y_arr)
        self._xiy2x = Transform2D.fit(xi_arr, y_arr, x_arr)
        self._xiy2lam = Transform2D.fit(xi_arr, y_arr, lam_arr)

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
        wave_min = fov.meta['wave_min'].value       # [um]
        wave_max = fov.meta['wave_max'].value       # [um]
        xi_min = fov.meta['xi_min'].value           # [arcsec]
        xi_max = fov.meta['xi_max'].value           # [arcsec]
        xlim_mm, ylim_mm = self.footprint(wave_min=wave_min, wave_max=wave_max,
                                          xi_min=xi_min, xi_max=xi_max)

        if xlim_mm is None:
            return None

        fov_header = fov.header
        det_header = fov.detector_header

        # WCSD from the FieldOfView - this is the full detector plane
        fpa_wcs = WCS(fov_header, key='D')
        naxis1, naxis2 = fov_header['NAXIS1'], fov_header['NAXIS2']
        pixsize = fov_header['CDELT1D'] * u.Unit(fov_header['CUNIT1D'])
        pixsize = pixsize.to(u.mm).value
        pixscale = fov_header['CDELT1'] * u.Unit(fov_header['CUNIT1'])
        pixscale = pixscale.to(u.arcsec).value

        fpa_wcsd = WCS(det_header, key='D')
        naxis1d, naxis2d = det_header['NAXIS1'], det_header['NAXIS2']
        xlim_px, ylim_px = fpa_wcsd.all_world2pix(xlim_mm, ylim_mm, 0)
        xmin = np.floor(xlim_px.min()).astype(int)
        xmax = np.ceil(xlim_px.max()).astype(int)
        ymin = np.floor(ylim_px.min()).astype(int)
        ymax = np.ceil(ylim_px.max()).astype(int)

        ## Check if spectral trace footprint is outside FoV
        if xmax < 0 or xmin > naxis1d or ymax < 0 or ymin > naxis2d:
            logging.warning("Spectral trace footprint is outside FoV")
            return None

        # Only work on parts within the FoV
        xmin = max(xmin, 0)
        xmax = min(xmax, naxis1d)
        ymin = max(ymin, 0)
        ymax = min(ymax, naxis2d)

        # Create header for the subimage - I think this only needs the DET one,
        # but we'll do both. The WCSs are initialised from the full fpa WCS and
        # then shifted accordingly.
        # sub_wcs = WCS(fov_header, key=" ")
        # sub_wcs.wcs.crpix -= np.array([xmin, ymin])
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

        # wavelength step per detector pixel at centre of slice
        # ..todo: I try the average here, but is that correct? The pixel
        #         edges may not correspond precisely to wave limits
        avg_dlam_per_pix = (wave_max - wave_min) / sub_naxis2
        try:
            xilam = XiLamImage(fov, avg_dlam_per_pix)
            self.xilam = xilam    # ..todo: remove
        except ValueError:
            print(" ---> ", self.meta['trace_id'], "gave ValueError")

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
            logging.warning("map_spectra_to_focal_plane:", np.sum(image < 0),
                            "negative pixels")

        image_hdu = fits.ImageHDU(header=img_header, data=image)
        return image_hdu


    def get_max_dispersion(self, **kwargs):
        '''Get the maximum dispersion in a spectral trace

        This is a wrapper for the function in `spectral_trace_utils`.
        '''
        params = {}
        params.update(self.meta)
        params["wave_min"] = self.wave_min
        params["wave_max"] = self.wave_max
        params.update(kwargs)

        disp, waverange = get_max_dispersion(self, **params)  # dwave is passed from kwargs
        self._disp = disp
        self._waverange = waverange

        return self._disp, self._waverange

    def get_pixel_wavelength_edges(self, pixel_size):
        """Returns the wavelengths at the edge of pixels along a trace"""
        self.meta["pixel_size"] = pixel_size
        if self._disp is None or self._waverange is None:
            self.get_max_dispersion()

        um_per_pix = pixel_size / self._disp  # wavelength range per pixel
        wbe = pixel_wavelength_edges(um_per_pix, self._waverange,
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
              ''.format(self.meta["trace_id"],
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
        cube_lam = wcs_lam.all_pix2world(np.arange(n_lam), 1)[0]
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

        self.image *= d_eta     # ph/s/um/arcsec2 --> ph/s/um/arcsec

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
        # ..todo: cubic spline introduces negative values, linear does not.
        #  Alternative might be to cubic-spline interpolate on sqrt(image),
        #  with subsequent squaring of the result. This would require
        #  wrapping RectBivariateSpline in a new (sub)class.
        spline_order = (1, 1)
        self.interp = RectBivariateSpline(self.xi, self.lam, self.image,
                                          kx=spline_order[0],
                                          ky=spline_order[1])
        # This is not executed. ..todo: define a switch?
        if False:
            fits.writeto("test_xilam.fits", data=self.image,
                         header=self.wcs.to_header(), overwrite=True)


class Transform2D():
    """
    2-dimensional polynomial transform

    The class is instantiated from a m x n matrix A that contains the
    coefficients of the polynomial. Along rows, the power of x increases;
    along columns, the power of y increases, such that A[j, i] is
    the coefficient of x^i y^j.

    Parameters
    ----------
    matrix : np.array
        matrix of polynomial coefficients

    ..todo: alternatively, the matrix can be created from a fit to data
    """

    def __init__(self, matrix):
        self.matrix = matrix
        self.ny, self.nx = matrix.shape

    def __call__(self, x, y, grid=False):
        """
        Apply the polynomial transform

        The transformation is a polynomial based on the simple
        monomials x^i y^j. When grid=True, the transform is applied to the grid
        formed by the tensor product of the vectors x and y. When grid=False,
        the vectors of x and y define the components of a number of points (in
        this case, x and y must be of the same length).

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
        x = np.array(x)
        y = np.array(y)
        orig_shape = x.shape

        if not grid and x.shape != y.shape:
            raise ValueError("x and y must have the same length when grid is False")

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
            if orig_shape == () or orig_shape is None:
                result = np.float32(result)
            else:
                result = result.reshape(orig_shape)

        return result

    @classmethod
    def fit(cls, xin, yin, xout, degree=4):
        """
        Determine polynomial fits
        """
        pinit = Polynomial2D(degree=degree)
        fitter = fitting.LinearLSQFitter()
        fit = fitter(pinit, xin, yin, xout)
        return Transform2D(fit2matrix(fit))

    def gradient(self):
        """Compute the gradient of a 2d polynomial transformation"""

        mat = self.matrix

        dmat_x = (mat * np.arange(self.nx))[:, 1:]
        dmat_y = (mat.T * np.arange(self.ny)).T[1:, :]

        return Transform2D(dmat_x), Transform2D(dmat_y)


def fit2matrix(fit):
    """
    Return coefficients from a polynomial fit as a matrix

    The Polynomial2D fits of degree n have coefficients for all i, j with i + j <= n.
    How would one rearrange those?
    """
    coeffs = dict(zip(fit.param_names, fit.parameters))
    deg = fit.degree
    mat = np.zeros((deg + 1 , deg + 1), dtype=np.float32)
    for i in range(deg + 1):
        for j in range(deg + 1):
            try:
                mat[j, i] = coeffs['c{}_{}'.format(i, j)]
            except KeyError:
                pass
    return mat

# ..todo: should the next three functions be combined and return a dictionary of fits?
def xilam2xy_fit(layout, params):
    """
    Determine polynomial fits of FPA position

    Fits are of degree 4 as a function of slit position and wavelength.
    """
    xi_arr = layout[params['s_colname']]
    lam_arr = layout[params['wave_colname']]
    x_arr = layout[params['x_colname']]
    y_arr = layout[params['y_colname']]

    ## Filter the lists: remove any points with x==0
    ## ..todo: this may not be necessary after sanitising the table
    #good = x != 0
    #xi = xi[good]
    #lam = lam[good]
    #x = x[good]
    #y = y[good]

    # compute the fits
    pinit_x = Polynomial2D(degree=4)
    pinit_y = Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xilam2x = fitter(pinit_x, xi_arr, lam_arr, x_arr)
    xilam2y = fitter(pinit_y, xi_arr, lam_arr, y_arr)

    return xilam2x, xilam2y

def xy2xilam_fit(layout, params):
    """
    Determine polynomial fits of wavelength/slit position

    Fits are of degree 4 as a function of focal plane position
    """

    xi_arr = layout[params['s_colname']]
    lam_arr = layout[params['wave_colname']]
    x_arr = layout[params['x_colname']]
    y_arr = layout[params['y_colname']]

    pinit_xi = Polynomial2D(degree=4)
    pinit_lam = Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xy2xi = fitter(pinit_xi, x_arr, y_arr, xi_arr)
    xy2lam = fitter(pinit_lam, x_arr, y_arr, lam_arr)

    return xy2xi, xy2lam


def _xiy2xlam_fit(layout, params):
    """Determine polynomial fits of wavelength/slit position

    Fits are of degree 4 as a function of focal plane position
    """

    # These are helper functions to allow fitting of left/right edges
    # for the purpose of checking whether a trace is on a chip or not.

    xi_arr = layout[params['s_colname']]
    lam_arr = layout[params['wave_colname']]
    x_arr = layout[params['x_colname']]
    y_arr = layout[params['y_colname']]

    pinit_x = Polynomial2D(degree=4)
    pinit_lam = Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xiy2x = fitter(pinit_x, xi_arr, y_arr, x_arr)
    xiy2lam = fitter(pinit_lam, xi_arr, y_arr, lam_arr)
    return xiy2x, xiy2lam


# ..todo: Check whether the following functions are actually used
def rolling_median(x, n):
    """ Calculates the rolling median of a sequence for +/- n entries """
    y = [np.median(x[max(0, i-n):min(len(x), i+n+1)]) for i in range(len(x))]
    return np.array(y)


def fill_zeros(x):
    """ Fills in zeros in a sequence with the previous non-zero number """
    for i in range(1, len(x)):
        if x[i] == 0:
            x[i] = x[i-1]
    return x

def get_max_dispersion(spt, wave_min, wave_max, dwave, **kwargs):
    """
    Finds the maximum distance [mm] per wavelength unit [um] along a trace

    The function looks for the minimum gradient of `spt.xy2lam` and returns
    the inverse and the corresponding wavelength.
    """
    gradient = deriv_polynomial2d(spt.xy2lam)
    dx_dwave = gradient[0](spt.table['x'], spt.table['y'])
    dy_dwave = gradient[1](spt.table['x'], spt.table['y'])
    absgrad2 = np.sqrt(dx_dwave**2 + dy_dwave**2)

    max_disp = 1. / np.min(absgrad2)
    wave_max = spt.table['wavelength'][np.argmin(absgrad2)]

    # ..todo: not sure why these have to be arrays?
    return np.array([max_disp]), np.array([wave_max])

def get_max_dispersion_old(trace_tbls, wave_min, wave_max, dwave,
                       **kwargs):
    """
    Finds the maximum distance [mm] per wavelength unit [um] along a trace

    Looks at all the trace lines (x, y) per wavelength across the slit for each
    trace, and for every trace projected onto the image plane.
    For each wavelength in the range [wave_min, wave_max] return the largest
    dx/dwave value based on all the trace projection lines.

    Parameters
    ----------
    trace_tbls : list of fits.BinTableHDU
        List of trace position tables. Units of table [um, mm, mm]
        Each table must have columns [wavelength, x0, y0, ..., xN, yN]
    wave_min, wave_max : float
        [um] minimum wavelength to look at
    dwave : float
        [um] wavelength step size

    kwargs
    ------
    x_colname, y_colname, wave_colname : str
        The name of each column for x, y, and wavelength on the image plane
        Default column names: ["x", "y", "wavelength"]
    col_number_start : int
        Default is 0. Start of the column numbering. I.e. x0, y0, s0 etc.
        If the columns start at x1, y1, etc; set ``col_number_start=1``

    Returns
    -------
    max_grad : array
        [mm/um] The maximum sensible gradient of all spectral trace projections
    waverange : array
        [um] The wavelengths corresponding to the gradients

    """

    params = {"x_colname": "x",
              "y_colname": "y",
              "wave_colname": "wavelength",
              "col_number_start": 0}
    params.update(kwargs)

    waverange = np.arange(wave_min, wave_max, dwave)
    dispersions = []
    for tbl in trace_tbls:
        if not isinstance(tbl, Table):
            tbl = Table(tbl)

        n = len([col for col in tbl.colnames if params["y_colname"] in col])
        k = params["col_number_start"]
        # .. todo: Can't use x1, etc. anymore, we have only one column x
        colnames = ["y"+str(ii) for ii in range(k, n+k)] + \
                   ["x"+str(ii) for ii in range(k, n+k)]
        for xcol in colnames:
            xpos = tbl[xcol]
            wave = tbl[params["wave_colname"]]
            if wave[0] > wave[1]:
                wave = wave[::-1]
                xpos = xpos[::-1]

            # disp is a new range [mm] derived from the trace coordinates (x, lam)
            # and the spectral resolution dwave
            mask = (waverange >= np.min(wave)) * (waverange <= np.max(wave))
            disp = np.zeros(len(waverange))
            disp[mask] = np.interp(waverange[mask], wave, xpos)
            disp /= dwave         # [mm/um] distance / wave_unit
            dispersions += [disp]

    # find the maximum dispersion for overlapping orders by using gradients
    # .. NOTE: careful of the np.abs(). Not sure if it should be here.
    grads = np.array([np.abs(np.gradient(disp)) for disp in dispersions])
    max_grad = fill_zeros(rolling_median(np.max(grads, axis=0), 15))

    # import matplotlib.pyplot as plt
    # for grad in grads:
    #     plt.plot(waverange, grad)
    # plt.scatter(waverange, max_grad)
    # plt.show()

    # max_grad is d_pos / d_wave : change in position [mm] per micron [um]
    return max_grad, waverange


def pixel_wavelength_edges(um_per_pix, waverange, wave_min, wave_max):
    """
    Get the wavelength bin edges for pixels under (a series) of spectral traces

    Returns the wavelength bin edges needed to properly project the spectrum
    according to the provided dispersion vector ``um_per_pix``

    Note: Units must be consistent, recommended [um]

    Parameters
    ----------
    um_per_pix : list, array
    waverange : list, array
    wave_min, wave_max : float

    Returns
    -------
    wave_bin_edges : array
        [um] The wavelength bin edges

    """

    wave_bin_edges = []
    wave = wave_min
    while wave < wave_max:
        wave_bin_edges += [wave]
        wave += np.interp(wave, waverange, um_per_pix)

    return np.array(wave_bin_edges)


def get_affine_parameters(coords):
    """
    Returns rotation and shear for each MTC point along a SpectralTrace

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
    shears = np.array([np.arctan2(dys[i], dxs[i]) for i in range(dxs.shape[0])])
    shears = np.array(list(shears.T) + [shears.T[-1]]).T
    shears = (np.average(shears, axis=0) * rad2deg) - (90 + rotations)

    return rotations, shears


def sanitize_table(tbl, invalid_value, wave_colname, x_colname, y_colname,
                   spline_order=4, ext_id=None):

    y_colnames = [col for col in tbl.colnames if y_colname in col]
    x_colnames = [col.replace(y_colname, x_colname) for col in y_colnames]

    for x_col, y_col in zip(x_colnames, y_colnames):
        wave = tbl[wave_colname].data
        x = tbl[x_col].data
        y = tbl[y_col].data

        valid = (x != invalid_value) * (y != invalid_value)
        invalid = np.invert(valid)
        if sum(invalid) == 0:
            continue

        if sum(valid) == 0:
            logging.warning("--- Extension {} ---"
                            "All points in {} or {} were invalid. \n"
                            "THESE COLUMNS HAVE BEEN REMOVED FROM THE TABLE \n"
                            "invalid_value = {} \n"
                            "wave = {} \nx = {} \ny = {}"
                            "".format(ext_id, x_col, y_col, invalid_value,
                                      wave, x, y))
            tbl.remove_columns([x_col, y_col])
            continue

        k = spline_order
        if wave[-1] > wave[0]:
            xnew = InterpolatedUnivariateSpline(wave[valid], x[valid], k=k)
            ynew = InterpolatedUnivariateSpline(wave[valid], y[valid], k=k)
        else:
            xnew = InterpolatedUnivariateSpline(wave[valid][::-1],
                                                x[valid][::-1], k=k)
            ynew = InterpolatedUnivariateSpline(wave[valid][::-1],
                                                y[valid][::-1], k=k)

        tbl[x_col][invalid] = xnew(wave[invalid])
        tbl[y_col][invalid] = ynew(wave[invalid])

    return tbl
