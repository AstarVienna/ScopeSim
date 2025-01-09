"""A class for the METIS WCU focal-plane mask"""

from pathlib import Path
import numpy as np
from astropy.io import fits
from ..data_container import DataContainer

class FPMask:
    """Focal-plane mask for the METIS WCU

    Parameters
    ----------
    See :class:`DataContainer` for input parameters

    """

    def __init__(self,
                 filename: Path | str | None = None,
                 **kwargs
                 ):
        self.filename = filename
        self.data_container = DataContainer(filename=filename, **kwargs)
        hdr = {"BG_SRC": True,
               "BG_SURF": "WCU focal plane mask",   # TODO more specific?
               "CTYPE1": "LINEAR",
               "CTYPE2": "LINEAR",
               "CRPIX1": 1024.5,
               "CRPIX2": 1024.5,
               "CRVAL1": 0.,
               "CRVAL2": 0.,
               "CUNIT1": "arcsec",
               "CUNIT2": "arcsec",
               "CDELT1": 0.00547,
               "CDELT2": 0.00547,
               "BUNIT": "PHOTLAM arcsec-2",
               "SOLIDANG": "arcsec-2"}

        self.hdu = self.make_hdu(header=hdr)


    def make_hdu(self, header) -> fits.ImageHDU:
        hdu = fits.ImageHDU()
        header = {"BG_SRC": True,    ## TODO this is a test
                  "BG_SURF": "WCU focal plane mask",   # TODO more specific?
                  "CTYPE1": "LINEAR",
                  "CTYPE2": "LINEAR",
                  "CRPIX1": 1024.,
                  "CRPIX2": 1024.,
                  "CRVAL1": 0.,
                  "CRVAL2": 0.,
                  "CUNIT1": "arcsec",
                  "CUNIT2": "arcsec",
                  "CDELT1": 0.00547,
                  "CDELT2": 0.00547,
                  "BUNIT": "PHOTLAM arcsec-2",
                  "SOLIDANG": "arcsec-2"}
        hdu.header.update(header)
        hdu.data = np.zeros((2047, 2047))  # TODO test - do we need to go back to 2048?
        tab = self.data_container.table
        xpix = (tab['x'] - header['CRVAL1']) / header['CDELT1'] + header['CRPIX1'] - 1
        ypix = (tab['y'] - header['CRVAL2']) / header['CDELT2'] + header['CRPIX2'] - 1
        holearea = (tab['diam']/2)**2 * np.pi
        hdu.data[ypix.astype(int), xpix.astype(int)] = holearea
        return hdu
