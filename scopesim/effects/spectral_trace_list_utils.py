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
    """Definition of one spectral trace

    A SpectralTrace describes the mapping of spectral slit coordinates
    to the focal plane. The class reads an order layout and fits several
    functions to describe the geometry of the trace.

    Slit coordinates are:
    - xi : spatial position along the slit [arcsec]
    - lam : Wavelength [um]
    Focal plane coordinates are:
    - x, y : [mm]
    """
    _class_params = {"x_colname": "x",
                     "y_colname": "y",
                     "s_colname": "s",
                     "wave_colname": "wavelength",
                     "dwave": 0.002,
                     "aperture_id": 0,
                     "image_plane_id": 0,
                     "extension_id": 2,
                     "spline_order": 4,
                     "pixel_size": None,
                     "description": "<no description>"}

    def __init__(self, trace_tbl, **kwargs):
        # Within scopesim, the actual parameter values are
        # passed as kwargs from the SpectralTraceList.
        # The values need to be here for stand-alone use.
        self.meta = {}
        self.meta.update(self._class_params)
        self.meta.update(kwargs)

        if isinstance(trace_tbl, (fits.BinTableHDU, fits.TableHDU)):
            self.table = Table.read(trace_tbl)
            self.meta["trace_id"] = trace_tbl.header.get('EXTNAME', "<unknown trace id>")
        elif isinstance(trace_tbl, Table):
            self.table = trace_tbl
        else:
            raise ValueError("trace_tbl must be one of (fits.BinTableHDU, "
                             "fits.TableHDU, astropy.Table): {}"
                             "".format(type(trace_tbl)))

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
        if self.meta["invalid_value"] is not None:
            self.table = sanitize_table(
                self.table,
                invalid_value=self.meta["invalid_value"],
                wave_colname=self.meta["wave_colname"],
                x_colname=self.meta["x_colname"],
                y_colname=self.meta["y_colname"],
                spline_order=self.meta["spline_order"],
                ext_id=self.meta["extension_id"])

        x_arr = self.table[self.meta['x_colname']]
        y_arr = self.table[self.meta['y_colname']]
        xi_arr = self.table[self.meta['s_colname']]
        lam_arr = self.table[self.meta['wave_colname']]

        wi0, wi1 = lam_arr.argmin(), lam_arr.argmax()
        x_disp_length = np.diff([x_arr[wi0], x_arr[wi1]])
        y_disp_length = np.diff([y_arr[wi0], y_arr[wi1]])
        self.dispersion_axis = "x" if x_disp_length > y_disp_length else "y"

        self.wave_min = quantify(np.min(lam_arr), u.um).value
        self.wave_max = quantify(np.max(lam_arr), u.um).value

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
        #print("xlim_mm:", xlim_mm, "   ylim_mm:", ylim_mm)
        if xlim_mm is None:
            print("xlim_mm is None")
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
        #print(fpa_wcsd)
        #print(xmin, xmax, ymin, ymax, " <<->> ", naxis1d, naxis2d)
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
        # ..todo: - currently using average dlam_per_pix. This should
        #           be okay if there is not strong anamorphism. Below, we
        #           compute an image of abs(dlam_per_pix) in the focal plane.
        #           XiLamImage would need that as an image of xi/lam, which should
        #           be possible but too much for the time being.
        #         - The dispersion direction is selected by the direction of the
        #           gradient of lam(x, y). This works if the lam-axis is well
        #           aligned with x or y. Needs to be tested for MICADO.


        # dlam_by_dx, dlam_by_dy = self.xy2lam.gradient()
        # if np.abs(dlam_by_dx(0, 0)) > np.abs(dlam_by_dy(0, 0)):
        if self.dispersion_axis == "x":
            avg_dlam_per_pix = (wave_max - wave_min) / sub_naxis1
        else:
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
            logging.warning(f"map_spectra_to_focal_plane: {np.sum(image < 0)} negative pixels")


        image_hdu = fits.ImageHDU(header=img_header, data=image)
        return image_hdu

    def footprint(self, wave_min=None, wave_max=None, xi_min=None, xi_max=None):
        """
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
        """
        #print(f"footprint: {wave_min}, {wave_max}, {xi_min}, {xi_max}")

        ## Define the wavelength range of the footprint. This is a compromise
        ## between the requested range (by method args) and the definition
        ## range of the spectral trace
        ## This is only relevant if the trace is given by a table of reference
        ## points. Otherwise (METIS LMS!) we assume that the range is valid.
        if ('wave_colname' in self.meta and
            self.meta['wave_colname'] in self.table.colnames):
            # Here, the parameters are obtained from a table of reference points
            wave_unit = self.table[self.meta['wave_colname']].unit
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

        return ([x_edge.min(), x_edge.max(), x_edge.max(), x_edge.min()],
                [y_edge.min(), y_edge.min(), y_edge.max(), y_edge.max()])

    def plot(self, wave_min=None, wave_max=None, c="r"):
        """Plot control points of the SpectralTrace"""

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
                plane_interp = RectBivariateSpline(cube_xi, cube_lam, plane,
                                                   kx=1, ky=1)
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


class Transform2D():
    """
    2-dimensional polynomial transform

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
        """Make sure `trafo` is a tuple"""
        if trafo is not None and not isinstance(trafo, tuple):
            trafo = (trafo, {})
        return trafo


    def __call__(self, x, y, grid=False, **kwargs):
        """
        Apply the polynomial transform

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
            raise ValueError("x and y must have the same length when grid is False")

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
            if orig_shape == () or orig_shape is None:
                result = np.float32(result)
            else:
                result = result.reshape(orig_shape)

        # Apply posttransform
        if self.posttransform is not None:
            result = self.posttransform[0](result, **self.posttransform[1])

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


# def sanitize_table(tbl, invalid_value, wave_colname, x_colname, y_colname,
#                    spline_order=4, ext_id=None):
#
#     y_colnames = [col for col in tbl.colnames if y_colname in col]
#     x_colnames = [col.replace(y_colname, x_colname) for col in y_colnames]
#
#     for x_col, y_col in zip(x_colnames, y_colnames):
#         wave = tbl[wave_colname].data
#         x = tbl[x_col].data
#         y = tbl[y_col].data
#
#         valid = (x != invalid_value) * (y != invalid_value)
#         invalid = np.invert(valid)
#         if sum(invalid) == 0:
#             continue
#
#         if sum(valid) == 0:
#             logging.warning("--- Extension {} ---"
#                             "All points in {} or {} were invalid. \n"
#                             "THESE COLUMNS HAVE BEEN REMOVED FROM THE TABLE \n"
#                             "invalid_value = {} \n"
#                             "wave = {} \nx = {} \ny = {}"
#                             "".format(ext_id, x_col, y_col, invalid_value,
#                                       wave, x, y))
#             tbl.remove_columns([x_col, y_col])
#             continue
#
#         k = spline_order
#         if wave[-1] > wave[0]:
#             xnew = InterpolatedUnivariateSpline(wave[valid], x[valid], k=k)
#             ynew = InterpolatedUnivariateSpline(wave[valid], y[valid], k=k)
#         else:
#             xnew = InterpolatedUnivariateSpline(wave[valid][::-1],
#                                                 x[valid][::-1], k=k)
#             ynew = InterpolatedUnivariateSpline(wave[valid][::-1],
#                                                 y[valid][::-1], k=k)
#
#         tbl[x_col][invalid] = xnew(wave[invalid])
#         tbl[y_col][invalid] = ynew(wave[invalid])
#
#     return tbl
