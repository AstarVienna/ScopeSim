import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy import units as u

from matplotlib import pyplot as plt

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
            self.table = Table(trace_tbl.data)
        elif isinstance(trace_tbl, Table):
            self.table = trace_tbl
        else:
            raise ValueError("trace_tbl must be one of (fits.BinTableHDU, "
                             "fits.TableHDU, astropy.Table): {}"
                             "".format(type(trace_tbl)))

        if self.meta["invalid_value"] is not None:
            self.table = spt_utils.sanitize_table(self.table,
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
        k, n = self.meta["col_number_start"], self.n_traces
        self.s_colnames = [self.meta["s_colname"]+str(i) for i in range(k, n+k)]
        self.x_colnames = [self.meta["x_colname"]+str(i) for i in range(k, n+k)]
        self.y_colnames = [self.meta["y_colname"]+str(i) for i in range(k, n+k)]

        self.wave_min = quantify(np.min(self.waves), u.um).value
        self.wave_max = quantify(np.max(self.waves), u.um).value

        self._disp = None           # [mm/um] spectral dispersion distance
        self._waverange = None
        self._wave_bin_edges = None
        self._wave_bin_centers = None
        self._curves = None

        self.xy2xi, self.xy2lam = spt_utils.xy2xilam_fit(self.table,
                                                         self.meta)
        self.xilam2x, self.xilam2y = spt_utils.xilam2xy_fit(self.table,
                                                            self.meta)
        self._xiy2x, self._xiy2lam = spt_utils._xiy2xlam_fit(self.table,
                                                             self.meta)

    def get_max_dispersion(self, **kwargs):
        params = {}
        params.update(self.meta)
        params["wave_min"] = self.wave_min
        params["wave_max"] = self.wave_max
        params.update(kwargs)

        disp, waverange = spt_utils.get_max_dispersion(trace_tbls=[self.table],
                                                       **params)  # dwave is passed from kwargs
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

    @property
    def footprint(self):
        '''Return corners of rectangle enclosing spectral trace'''
        xval = self.table[self.meta['x_colname']]
        yval = self.table[self.meta['y_colname']]
        xlim = [np.min(xval), np.max(xval), np.max(xval), np.min(xval)]
        ylim = [np.min(yval), np.min(yval), np.max(yval), np.max(yval)]
        return xlim, ylim

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
        xlim, ylim  = self.footprint
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
