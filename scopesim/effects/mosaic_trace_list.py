#  -*- coding: utf-8 -*-
"""SpectralTraceList and SpectralTrace for MOSAIC"""
from tqdm.auto import tqdm

import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.wcs import WCS
from astropy.modeling import fitting
from astropy.modeling.models import Polynomial1D

from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace

from ..utils import get_logger, quantify, power_vector
from ..optics.fov import FieldOfView
from ..optics.fov_volume_list import FovVolumeList

logger = get_logger(__name__)

class MosaicSpectralTraceList(SpectralTraceList):
    """SpectralTraceList for MOSAIC"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.aplist = self._file["Aperture List"].data
        self.view = np.array([(self.aplist["right"].max() -
                               self.aplist["left"].min()),
                              (self.aplist["top"].max() -
                               self.aplist["bottom"].min())])

    def apply_to(self, obj, **kwargs):
        """See parent docstring."""
        ### This is copied from MetisSpectralTraceList, make less redundant?
        if isinstance(obj, FovVolumeList):
            logger.debug("Executing %s, FoV setup", self.meta['name'])
            # Create a single volume that covers the aperture and
            # the maximum wavelength range of the grating
            volumes = [spectral_trace.fov_grid()
                       for spectral_trace in self.spectral_traces.values()]
            wave_min = min(vol["wave_min"] for vol in volumes)
            wave_max = max(vol["wave_max"] for vol in volumes)
            extracted_vols = obj.extract(axes=["wave"],
                                         edges=([[wave_min, wave_max]]))
            obj.volumes = extracted_vols
            print("Mosaic Spectral Trace List:", obj)

        if isinstance(obj, FieldOfView):
            # Application to field of view
            logger.debug("Executing %s, FoV", self.meta['name'])
            if obj.hdu is not None and obj.hdu.header["NAXIS"] == 3:
                obj.cube = obj.hdu
            elif obj.hdu is None and obj.cube is None:
                print("MosSpTrL: Making cube")
                obj.cube = obj.make_cube_hdu()

            fovcube = obj.cube.data
            n_z, n_y, n_x = fovcube.shape
            fovwcs = WCS(obj.cube.header)
            # Make this linear to avoid jump at RA 0 deg
            fovwcs.wcs.ctype = ["LINEAR", "LINEAR", fovwcs.wcs.ctype[2]]
            fovwcs_spat = fovwcs.sub(2)

            ## This is the place where we need to look at the apertures
            ## - collapse each aperture to 1D spectrum by integrating spatially
            ## - map each 1D spectrum to detector/fov
            fovimage = np.zeros((obj.detector_header["NAXIS2"],
                                 obj.detector_header["NAXIS1"]),
                                dtype=np.float32)

            for sptid, spt in tqdm(self.spectral_traces.items(),
                                   desc="Fiber traces", position=2):
                nx_slice = (self.aplist[sptid]["right"] - self.aplist[sptid]["left"])/pixscale
                nx_slice = (self.aplist[sptid]["top"] - self.aplist[sptid]["bottom"])/pixscale

                ymin = spt.meta["fov"]["y_min"]
                ymax = spt.meta["fov"]["y_max"]

                slicewcs = fovwcs.deepcopy()
                slicewcs.wcs.ctype = ["LINEAR", "LINEAR",
                                      slicewcs.wcs.ctype[2]]
                slicewcs.wcs.crpix[0] = (nx_slice + 1) / 2
                slicewcs.wcs.crpix[1] = (ny_slice + 1) / 2
                slicewcs.wcs.crval[0] = (xmin + xmax) / 2 / 3600
                slicewcs.wcs.crval[1] = (ymin + ymax) / 2 / 3600
                slicewcs.wcs.cdelt[0] = (ymax - ymin) / ny_slice / 3600
                slicewcs.wcs.cdelt[1] = (ymax - ymin) / ny_slice / 3600
                slicewcs_spat = slicewcs.sub(2)

                # World coordinates for the slice
                xworld, yworld = slicewcs_spat.app_pix2world(xslice, yslice, 0)
                # FOV pixel coordinates for the slice
                xfov, yfov = fovwcs_spat.all_world2pix(xworld, yworld, 0)

                for islice in range(n_z):
                    # this should not be necessary
                    ifov = RectBivariateSpline(np.arange(n_y),
                                               np.arange(n_x),
                                               fovcube[islice], kx=1, ky=1)
                    fovcube[islice, yfov, xfov].sum()
                # build slicefov = FieldOfView1D(...)
                # fovimage[] += slicefov.hdu.data (better integer row = slicefov...)

        return obj



    def make_spectral_traces(self):
        """Return a dictionary of spectral traces read in from a file."""
        self.ext_data = self._file[0].header["EDATA"]
        self.ext_cat = self._file[0].header["ECAT"]
        self.catalog = Table(self._file[self.ext_cat].data)
        spec_traces = {}
        for row in self.catalog:
            if row["image_plane_id"] == -1:
                continue
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            hdu = self._file[row["extension_id"]]
            spec_traces[row["description"]] = MosaicSpectralTrace(hdu, **params)

        self.spectral_traces = spec_traces


class MosaicSpectralTrace(SpectralTrace):

    def __init__(self, trace_tbl, **kwargs):
        super().__init__(trace_tbl, **kwargs)



    def compute_interpolation_functions(self):
        x_arr = self.table[self.meta["x_colname"]]
        y_arr = self.table[self.meta["y_colname"]]
        xi_arr = self.table[self.meta["s_colname"]]
        lam_arr = self.table[self.meta["wave_colname"]]

        self.wave_min = quantify(np.min(lam_arr), u.um).value
        self.wave_max = quantify(np.max(lam_arr), u.um).value



class Transform1D():
    """
    1-dimensional polynomial transform.
    """

    def __init__(self, coeffs, pretransform=None,
                 posttransform=None):
        self.coeffs = np.asarray(coeffs)
        self.nx = self.coeffs.shape[0]
        self.pretransform = self._repackage(pretransform)
        self.posttransform = self._repackage(posttransform)

    def _repackage(self, trafo):
        """Make sure `trafo` is a tuple."""
        if trafo is not None and not isinstance(trafo, tuple):
            trafo = (trafo, {})
        return trafo

    def __call__(self, x, **kwargs):
        """
        Apply the polynomial transform.

        The transformation is a polynomial based on the simple monomials x^i.
        """

        if "pretransform" in kwargs:
            self.pretransform = self._repackage(kwargs["pretransform"])
        if "postransform" in kwargs:
            self.posttransform = self._repackage(kwargs["posttransform"])

        x = np.array(x)

        # Apply pre transform
        if self.pretransform is not None:
            x = self.pretransform[0](x, **self.pretransform[1])

        xvec = power_vector(x, self.nx - 1)

        result = self.coeffs @ xvec

        # Apply posttransform
        if self.posttransform is not None:
            result = self.posttransform[0](result, **self.posttransform[1])

        return result

    @classmethod
    def fit(cls, xin, xout, degree=4):
        """Determine polynomial fit"""
        pinit = Polynomial1D(degree=degree)
        fitter = fitting.LinearLSQFitter()
        fit = fitter(pinit, xin, xout)
        return Transform1D(fit.parameters)

    def gradient(self):
        """Compute the gradient of a 1d polynomial transformation"""
        coeffs = self.coeffs

        dcoeffs = (coeffs * np.arange(self.nx))[1:]
        return Transform1D(dcoeffs)
