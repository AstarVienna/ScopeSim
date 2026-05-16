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
from .spectral_trace_list_utils import make_image_interpolations
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
    }


    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # field of view of the instrument

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
            hdu = self._file[row["extension_id"]]
            for apid in self.slicelist["id"]:
                specid = f"{row['description']}_{apid:02d}"
                params["aperture_id"] = apid
                spec_traces[specid] = SpectralTrace(hdu, **params)

        self.spectral_traces = spec_traces
