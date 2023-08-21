"""
Spectral grating efficiencies
"""
import logging
import numpy as np
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table

from .effects import Effect
from .ter_curves import TERCurve
from .ter_curves_utils import apply_throughput_to_cube
from ..utils import find_file
from ..base_classes import FieldOfViewBase, FOVSetupBase

class SpectralEfficiency(Effect):
    """
    Applies the grating efficiency (blaze function) for a SpectralTraceList

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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        params = {"z_order": [630]}
        self.meta.update(params)

        self.efficiencies = self.get_efficiencies()

    def get_efficiencies(self):
        """Reads effciencies from file, returns a dictionary"""
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
            effic_curve = TERCurve(wavelength=wavelength,
                                   transmission=efficiency,
                                   **params)
            efficiencies[name] = effic_curve

        hdul.close()
        return efficiencies

    def apply_to(self, obj, **kwargs):
        """
        Interface between FieldOfView and SpectralEfficiency
        """
        trace_id = obj.meta['trace_id']
        try:
            effic = self.efficiencies[trace_id]
        except KeyError:
            logging.warning("No grating efficiency for trace %s" % trace_id)
            return obj

        wcs = WCS(obj.hdu.header).spectral
        wave_cube = wcs.all_pix2world(np.arange(obj.hdu.data.shape[0]), 0)[0]
        wave_cube = (wave_cube * u.Unit(wcs.wcs.cunit[0])).to(u.AA)
        obj.hdu = apply_throughput_to_cube(obj.hdu, effic.throughput)
        return obj

    def plot(self):
        """Plot the grating efficiencies"""
        for name, effic in self.efficiencies.items():
            wave = effic.throughput.waveset
            plt.plot(wave.to(u.um), effic.throughput(wave), label=name)

        plt.xlabel("Wavelength [um]")
        plt.ylabel("Grating efficiency")
        plt.title(f"Grating efficiencies from {self.filename}")
        plt.legend()
        plt.show()
