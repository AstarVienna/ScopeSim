# -*- coding: utf-8 -*-
"""Spectral grating efficiencies."""

from typing import ClassVar

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table

from .effects import Effect
from .ter_curves import TERCurve
from .ter_curves_utils import apply_throughput_to_cube
from ..utils import figure_factory, get_logger


logger = get_logger(__name__)


class SpectralEfficiency(Effect):
    """
    Applies the grating efficiency (blaze function) for a SpectralTraceList.

    Input Data Format
    -----------------
    The efficiency curves are taken from a fits file `filename`with a
    structure similar to the trace definition file (see `SpectralTraceList`).
    The required extensions are:
    - 0 : PrimaryHDU [header]
    - 1 : BinTableHDU or TableHDU[header, data] : Overview table of all traces
    - 2..N : BinTableHDU or TableHDU : Efficiency curves, one per trace. The
             tables must have the two columns `wavelength` and `efficiency`

    Note that there must be one extension for each trace defined in the
    `SpectralTraceList`. Extensions for other traces are ignored.

    EXT 0 : PrimaryHDU
    ++++++++++++++++++
    Required header keywords:

    - ECAT : int : Extension number of overview table, normally 1
    - EDATA : int : Extension number of first Trace table, normally 2

    No data is required in this extension

    EXT 1 : (Bin)TableHDU : Overview of traces
    ++++++++++++++++++++++++++++++++++++++++++
    No special header keywords are required in this extension.

    Required Table columns:
    - description : str : identifier for each trace
    - extension_id : int : which extension is each trace in

    EXT 2 : (Bin)TableHDU : Efficiencies for individual traces
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Required header keywords:
    - EXTNAME : must be identical to the `description` in EXT 1

    Required Table columns:
    - wavelength : float : [um]
    - efficiency : float : number [0..1]

    """

    z_order: ClassVar[tuple[int, ...]] = (630,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        self.efficiencies = self.get_efficiencies()

    def get_efficiencies(self):
        """Read effciencies from file, returns a dictionary."""
        hdul = self._file
        self.ext_data = hdul[0].header["EDATA"]
        self.ext_cat = hdul[0].header["ECAT"]
        self.catalog = Table(hdul[self.ext_cat].data)

        efficiencies = {}
        for row in self.catalog:
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            hdu = self._file[row["extension_id"]]
            name = hdu.header['EXTNAME']

            tbl = Table.read(hdu)
            wavelength = tbl['wavelength'].quantity
            efficiency = tbl['efficiency'].value
            params.pop("filename", None)  # don't pass filename to TERCurve!
            effic_curve = TERCurve(array_dict={"wavelength":wavelength,
                                   "transmission":efficiency},
                                   **params)
            efficiencies[name] = effic_curve

        hdul.close()
        return efficiencies

    def apply_to(self, obj, **kwargs):
        """Interface between FieldOfView and SpectralEfficiency."""
        trace_id = obj.trace_id
        try:
            effic = self.efficiencies[trace_id]
        except KeyError:
            logger.warning("No grating efficiency for trace %s", trace_id)
            return obj

        swcs = WCS(obj.hdu.header).spectral
        with u.set_enabled_equivalencies(u.spectral()):
            wave = swcs.pixel_to_world(np.arange(swcs.pixel_shape[0])) << u.um
        obj.hdu = apply_throughput_to_cube(obj.hdu, effic.throughput, wave)
        return obj

    def plot(self):
        """Plot the grating efficiencies."""
        fig, axes = figure_factory()
        for name, effic in self.efficiencies.items():
            wave = effic.throughput.waveset
            axes.plot(wave.to(u.um), effic.throughput(wave), label=name)

        axes.set_xlabel("Wavelength [um]")
        axes.set_ylabel("Grating efficiency")
        axes.set_title(f"Grating efficiencies {self.display_name}")
        axes.legend()

        return fig
