#  -*- coding: utf-8 -*-
"""SpectralTraceList and SpectralTrace for MOSAIC"""
import numpy as np
from astropy.table import Table
from astropy import units as u

from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace

from ..utils import get_logger, quantify

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
        print(self.table)
        x_arr = self.table[self.meta["x_colname"]]
        y_arr = self.table[self.meta["y_colname"]]
        xi_arr = self.table[self.meta["s_colname"]]
        lam_arr = self.table[self.meta["wave_colname"]]

        self.wave_min = quantify(np.min(lam_arr), u.um).value
        self.wave_max = quantify(np.max(lam_arr), u.um).value
