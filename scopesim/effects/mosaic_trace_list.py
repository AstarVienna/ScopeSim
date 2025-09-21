#  -*- coding: utf-8 -*-
"""SpectralTraceList and SpectralTrace for MOSAIC."""

from typing import ClassVar

from tqdm.auto import tqdm
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling import fitting
from astropy.modeling.models import Polynomial1D
from synphot import SourceSpectrum, Empirical1D

from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace
from ..utils import get_logger, quantify, power_vector
from ..optics.fov import FieldOfView
from ..optics.fov_volume_list import FovVolumeList
from ..detector import Detector

logger = get_logger(__name__)


class MosaicSpectralTraceList(SpectralTraceList):
    """SpectralTraceList for MOSAIC.

    .. versionadded:: 0.11.0

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.aplist = self._file["Aperture List"].data
        # TODO: check units or normalise to arcsec
        self.view = np.array([(self.aplist["right"].max() -
                               self.aplist["left"].min()),
                              (self.aplist["top"].max() -
                               self.aplist["bottom"].min())])

    def apply_to(self, obj, **kwargs):
        """See parent docstring."""
        # This is copied from MetisSpectralTraceList, make less redundant?
        if isinstance(obj, FovVolumeList):
            logger.debug("Executing %s, FoV setup", self.meta["name"])
            # Create a single volume that covers the aperture and
            # the maximum wavelength range of the grating
            volumes = [spectral_trace.fov_grid()
                       for spectral_trace in self.spectral_traces.values()]
            wave_min = min(vol["wave_min"] for vol in volumes)
            wave_max = max(vol["wave_max"] for vol in volumes)
            extracted_vols = obj.extract(axes=["wave"],
                                         edges=([[wave_min, wave_max]]))
            obj.volumes = extracted_vols

        if isinstance(obj, FieldOfView):
            # Application to field of view
            logger.debug("Executing %s, FoV", self.meta["name"])
            if obj.hdu is not None and obj.hdu.header["NAXIS"] == 3:
                obj.cube = obj.hdu
            elif obj.hdu is None and obj.cube is None:
                obj.cube = obj.make_cube_hdu()

            fovcube = obj.cube.data
            n_z = fovcube.shape[0]
            fovwcs = WCS(obj.cube.header)
            # Make this linear to avoid jump at RA 0 deg
            fovwcs.wcs.ctype = ["LINEAR", "LINEAR", fovwcs.wcs.ctype[2]]
            fovwcs_spat = fovwcs.sub(2)
            fovwcs_spec = fovwcs.spectral
            fovlam = fovwcs_spec.all_pix2world(np.arange(n_z), 0)[0]
            fovlam <<= u.Unit(fovwcs_spec.wcs.cunit[0])

            det_header = obj.detector_header
            detwcs = WCS(det_header, key="D")
            naxis1d, naxis2d = det_header["NAXIS1"], det_header["NAXIS2"]

            # This is the place where we need to look at the apertures
            # - collapse each aperture to 1D spectrum by integrating spatially
            # - map each 1D spectrum to detector/fov

            image = np.zeros((naxis2d, naxis1d), dtype=np.float32)

            for sptid, spt in tqdm(self.spectral_traces.items(),
                                   desc="Fiber traces", position=2):
                theap = self.aplist[self.aplist["id"] == sptid]

                # solid angle in arcsec**2
                solid_angle = ((theap["right"] - theap["left"]) *
                               (theap["top"] - theap["bottom"]))

                # apertures are defined in arcsec. fovwcs is in degrees
                xmin, xmax, ymin, ymax = (theap["left"]/3600, theap["right"]/3600,
                                          theap["bottom"]/3600, theap["top"]/3600)

                imin = max(0, int(np.round(fovwcs_spat.all_world2pix(xmin, 0, 0)[0][0])))
                imax = int(np.round(fovwcs_spat.all_world2pix(xmax, 0, 0)[0][0]))
                jmin = max(0, int(np.round(fovwcs_spat.all_world2pix(0, ymin, 0)[1][0])))
                jmax = int(np.round(fovwcs_spat.all_world2pix(0, ymax, 0)[1][0]))

                # Average over the spatial dimensions of the aperture (still per arcsec2)
                fovflux = fovcube[:, jmin:jmax, imin:imax].mean(axis=(1,2)) * solid_angle
                spec = SourceSpectrum(Empirical1D, points=fovlam.to(u.um),
                                      lookup_table=fovflux)

                # Need to interpolate this to the output wavelength grid
                detlam = spt.x2lam(detwcs.all_pix2world(np.arange(naxis1d), 0, 0)[0])
                detlam <<= u.um
                yvals = spt.lam2y(detlam.value)
                jfib = detwcs.all_world2pix(0, yvals.mean(), 0)[1].astype(int)
                logger.debug("Flux from %s: %.4g", spt.trace_id, spec(detlam).value.sum())

                detdisp = np.diff(detlam, prepend=detlam[0])
                image[jfib,] += (spec(detlam) * detdisp).value

            image_hdr = detwcs.to_header()
            image_hdr["BUNIT"] = "ph s-1"
            image_hdr.extend(det_header)
            obj.hdu = fits.ImageHDU(data=image, header=image_hdr)
        return obj

    def make_spectral_traces(self):
        """Return a dictionary of spectral traces read in from a file."""
        self.ext_data = self._file[0].header["EDATA"]
        self.ext_cat = self._file[0].header["ECAT"]
        self.catalog = Table(self._file[self.ext_cat].data)
        spec_traces = {}
        for row in self.catalog:
            # image_plane_id = -1 marks rows that should not be read,
            # e.g. the aperture list. Although not necessary if the catalogue
            # is formatted in a way that only traces are listed, this provides
            # a possibility to "mask" traces.
            if row["image_plane_id"] == -1:
                continue
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            hdu = self._file[row["extension_id"]]
            spec_traces[row["description"]] = MosaicSpectralTrace(hdu, **params)

        self.spectral_traces = spec_traces


class MosaicSpectralTrace(SpectralTrace):
    """A single spectral trace for MOSAIC.

    .. versionadded:: 0.11.0
    """

    def compute_interpolation_functions(self):
        x_arr = self.table[self.meta["x_colname"]]
        y_arr = self.table[self.meta["y_colname"]]
        # xi_arr = self.table[self.meta["s_colname"]]
        lam_arr = self.table[self.meta["wave_colname"]]

        self.wave_min = quantify(np.min(lam_arr), u.um).value
        self.wave_max = quantify(np.max(lam_arr), u.um).value

        self.lam2x = Transform1D.fit(lam_arr, x_arr, degree=2)
        self.x2lam = Transform1D.fit(x_arr, lam_arr, degree=2)
        self.lam2y = Transform1D.fit(lam_arr, y_arr, degree=2)


class Transform1D:
    """1-dimensional polynomial transform.

    .. versionadded:: 0.11.0
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
    def fit(cls, xin, xout, degree: int = 4):
        """Determine polynomial fit."""
        pinit = Polynomial1D(degree=degree)
        fitter = fitting.LinearLSQFitter()
        fit = fitter(pinit, xin, xout)
        return Transform1D(fit.parameters)

    def gradient(self):
        """Compute the gradient of a 1d polynomial transformation."""
        coeffs = self.coeffs

        dcoeffs = (coeffs * np.arange(self.nx))[1:]
        return Transform1D(dcoeffs)


class MosaicCollapseSpectralTraces(MosaicSpectralTraceList):
    """Collapse SpectralTraces to 1D spectrum.

    .. versionadded:: 0.11.0
    """

    required_keys = {"filename"}
    z_order: ClassVar[tuple[int, ...]] = (899,)

    def apply_to(self, det, **kwargs):
        """Apply to detector readout."""
        if not isinstance(det, Detector):
            return det

        image = det._hdu.data
        detwcs = WCS(det._hdu.header, key="D")
        spec = np.zeros(image.shape[1], dtype=np.float32)
        for sptid, spt in tqdm(self.spectral_traces.items(),
                               desc="Fiber traces", position=2):
            y_mm = spt.table["y"][0]
            jfib = int(detwcs.all_world2pix(0, y_mm, 0)[1])
            spec += image[jfib,]

        x_mm = detwcs.all_pix2world(np.arange(image.shape[1]), 1, 0)[0]
        lam = spt.x2lam(x_mm)
        det._hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="wavelength", format="D", array=lam, unit="um"),
            fits.Column(name="spectrum", format="D", array=spec, unit="ADU"),
        ])

        return det
