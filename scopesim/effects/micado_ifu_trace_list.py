# -*- coding: utf-8 -*-
"""SpectralTraceList and SpectralTrace for the MICADO IFU."""

import copy
import warnings

from tqdm.auto import tqdm
import numpy as np
from scipy.interpolate import RectBivariateSpline

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u

from ..utils import from_currsys, find_file, quantify, get_logger
from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace
from .spectral_trace_list_utils import Transform2D
from .spectral_trace_list_utils import make_image_interpolations, det_offset
from .apertures import ApertureMask
from .ter_curves import TERCurve
from ..optics.fov import FieldOfView, FieldOfView3D
from ..optics.fov_volume_list import FovVolumeList

logger = get_logger(__name__)


class MicadoIFUSpectralTraceList(SpectralTraceList):
    """SpectralTraceList for the MICADO IFU."""

    _class_params = {
        "naxis1": 112,
        "nslice": 32,
        "slicewidth": 0.012, # arcsec
        "wave_min": 1.5,
        "wave_max": 2.5,   # todo: do not hardcode
    }


    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def apply_to_fovvolumelist(self, obj):
        """
        Setup

        Create a single volume that covers the aperture and the maximum
        wavelength range of the MICADO IFU.
        """
        logger.debug("Executing %s: Setup", self.meta['name'])
        volumes = [spectral_trace.fov_grid()
                   for spectral_trace in self.spectral_traces.values()]

        wave_min = min(vol["wave_min"] for vol in volumes)
        wave_max = max(vol["wave_max"] for vol in volumes)
        extracted_vols = obj.extract(axes=["wave"],
                                     edges=([[wave_min, wave_max]]))
        obj.volumes = extracted_vols
        self.wave_min = wave_min
        self.wave_max = wave_max
        return obj

    def apply_to_fov(self, obj):
        """Apply to FieldOfView"""

        logger.debug("Executing %s, FoV", self.meta['name'])
        if obj.hdu is not None and obj.hdu.header["NAXIS"] == 3:
            obj.cube = obj.hdu
        elif obj.hdu is None and obj.cube is None:
            obj.cube = obj.make_hdu()

        fovcube = obj.cube.data
        n_z, n_y, n_x = fovcube.shape
        fovwcs = WCS(obj.cube.header)
        # Make this linear to avoid jump at RA 0 deg
        fovwcs.wcs.ctype = ["LINEAR", "LINEAR", fovwcs.wcs.ctype[2]]
        fovwcs_spat = fovwcs.sub(2)
        ny_slice = self.meta["slice_samples"]

        # Spatial pixel coordinates
        xslice, yslice = np.meshgrid(np.arange(n_x),
                                     np.arange(ny_slice))

        fovimage = np.zeros((obj.detector_header["NAXIS2"],
                             obj.detector_header["NAXIS1"]),
                            dtype=np.float32)

        for sptid, spt in tqdm(self.spectral_traces.items(),
                               desc=" Spectral Traces", position=2):
            ymin = spt.meta["fov"]["y_min"]
            ymax = spt.meta["fov"]["y_max"]


            slicewcs = fovwcs.deepcopy()

            slicewcs.wcs.ctype = ["LINEAR", "LINEAR",
                                  slicewcs.wcs.ctype[2]]
            slicewcs.wcs.crpix[1] = (ny_slice + 1) / 2
            slicewcs.wcs.crval[1] = (ymin + ymax) / 2 / 3600
            slicewcs.wcs.cdelt[1] = (ymax - ymin) / ny_slice / 3600
            slicewcs_spat = slicewcs.sub(2)

            # World coordinates for the slice
            xworld, yworld = slicewcs_spat.all_pix2world(xslice,
                                                         yslice, 0)
            # FOV pixel coordinates for the slice
            xfov, yfov = fovwcs_spat.all_world2pix(xworld, yworld, 0)

            slicecube = np.zeros((n_z, ny_slice, n_x))
            for islice in range(n_z):
                ifov = RectBivariateSpline(np.arange(n_y),
                                           np.arange(n_x),
                                           fovcube[islice], kx=1, ky=1)
                slicecube[islice] = ifov(yfov, xfov, grid=False)

            slicefov = FieldOfView3D(obj.header,
                                     [obj.meta["wave_min"],
                                      obj.meta["wave_max"]])
            slicefov.detector_header = obj.detector_header
            slicefov.meta["xi_min"] = spt.meta["fov"]["xi_min"]
            slicefov.meta["xi_max"] = spt.meta["fov"]["xi_max"]
            slicefov.meta["trace_id"] = sptid
            slicefov.cube = fits.ImageHDU(header=slicewcs.to_header(),
                                          data=slicecube)
            # slicefov.cube.writeto(f"slicefov_{sptid}.fits", overwrite=True)
            slicefov.hdu = spt.map_spectra_to_focal_plane(slicefov)
            if slicefov.hdu is not None:
                sxmin = slicefov.hdu.header["XMIN"]
                sxmax = slicefov.hdu.header["XMAX"]
                symin = slicefov.hdu.header["YMIN"]
                symax = slicefov.hdu.header["YMAX"]
                fovimage[symin:symax, sxmin:sxmax] += slicefov.hdu.data

        obj.hdu = fits.ImageHDU(data=fovimage, header=obj.detector_header)
        return obj




    def make_spectral_traces(self):
        """Make a spectral trace for each combination of order and aperture"""
        self.ext_data = self._file[0].header["EDATA"]
        self.ext_cat = self._file[0].header["ECAT"]
        self.catalog = Table(self._file[self.ext_cat].data)
        self.slicelist = self._file["Aperture List"].data
        spec_traces = {}
        for row in self.catalog:
            if row["image_plane_id"] == -99:
                continue
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            trace_hdu = self._file[row["extension_id"]]

            for spslice, apid in enumerate(self.slicelist["id"]):
                specid = f"{row['description']}_{apid:02d}"
                thistrace = MicadoIFUSpectralTrace(trace_hdu, self.slicelist, spslice,
                                                   **params)
                thistrace.meta["trace_id"] = specid
                spec_traces[specid] = thistrace

        self.spectral_traces = spec_traces



class MicadoIFUSpectralTrace(SpectralTrace):
    """SpectralTrace for the MICADO IFU"""

    _class_params = {
        "naxis1": 112,
        "nslice": 32,
        "slicewidth": 0.012, # arcsec
        "pixsize": 0.015,    # mm
    }


    def __init__(self, trace_tbl, aplist, spslice, cmds=None, **kwargs):
        super().__init__(trace_tbl, **kwargs)

        self._set_dispersion(self.wave_min, self.wave_max)
        # Provisional
        self.aplist = aplist
        self.meta["slice"] = spslice
        self.meta["fov"] = self.fov_grid()

        lam = np.linspace(self.wave_min, self.wave_max, 1001) * u.um
        off_y = aplist['offset'][spslice]
        dlam = np.median(self.dlam_per_pix(lam)) / self.meta["pixsize"] * off_y
        self.xy2lam.posttransform = (det_offset, {"offset": dlam})
        self.xilam2x.pretransform_y = (det_offset, {"offset": dlam})
        self.xilam2y.pretransform_y = (det_offset, {"offset": dlam})

    def fov_grid(self):
        """
        Provide information on the source space volume required by the effect.
        """
        aperture = self.aplist[self.meta["slice"]]
        x_min = aperture["left"]
        x_max = aperture["right"]
        y_min = aperture["bottom"]
        y_max = aperture["top"]
        xi_min = aperture["xi1"]
        xi_max = aperture["xi2"]

        # ..todo: just a hack - xi and x are the same except xi is a quantity
        xi_min = quantify(xi_min, u.arcsec)
        xi_max = quantify(xi_max, u.arcsec)

        return {"x_min": x_min, "x_max": x_max,
                "y_min": y_min, "y_max": y_max,
                "xi_min": xi_min, "xi_max": xi_max,
                "wave_min": self.meta["wave_min"],
                "wave_max": self.meta["wave_max"],
                "trace_id": self.trace_id}
