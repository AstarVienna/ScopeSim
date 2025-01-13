"""A class for the METIS WCU focal-plane mask"""

from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from ..data_container import DataContainer
from ...utils import find_file, from_currsys, get_logger

logger = get_logger(__name__)

class FPMask:
    """Focal-plane mask for the METIS WCU

    Parameters
    ----------
    See :class:`DataContainer` for input parameters

    """

    def __init__(self,
                 maskname: Path | str | None = None,
                 fpmask_filename_format: str | None = None,
                 **kwargs
                 ):
        logger.debug("Initialising FPMask with ", maskname)
        # Try to find the file as a path
        if find_file(maskname, silent=True) is None:
            file_format = from_currsys(fpmask_filename_format)
            self.filename = file_format.format(maskname)
        else:
            self.filename = maskname

        self.data_container = DataContainer(filename=self.filename, **kwargs)
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

        self.pixarea = (hdr['CDELT1'] * u.Unit(hdr['CUNIT1'])
                        * hdr['CDELT2'] * u.Unit(hdr['CUNIT2']))

        self.holehdu = self.make_holehdu(header=hdr)
        self.opaquehdu = self.make_opaquehdu(header=hdr)


    def make_holehdu(self, header) -> fits.ImageHDU:
        """Create an hdu for the holes in fpmask

        The holes are assumed to be unresolved. They therefore cover one pixel and have
        a value corresponding to the actual solid angle covered by the hole.
        """

        hdu = fits.ImageHDU()
        hdu.header.update(header)
        hdu.data = np.zeros((2047, 2047))  # TODO test - do we need to go back to 2048?
        tab = self.data_container.table
        xpix = (tab['x'] - header['CRVAL1']) / header['CDELT1'] + header['CRPIX1'] - 1
        ypix = (tab['y'] - header['CRVAL2']) / header['CDELT2'] + header['CRPIX2'] - 1
        holearea = (tab['diam']/2)**2 * np.pi
        hdu.data[ypix.astype(int), xpix.astype(int)] = holearea
        return hdu

    def make_opaquehdu(self, header) -> fits.ImageHDU:
        """Create an hdu for the opaque area of the mask

        When spectrum is intensity (../arcsec2), the image must contain the pixel
        area (arcsec2). For a true backgroud field, this is taken care of by
        fov._calc_area_factor, but this is not applied to image data, so we
        have to do it here.

        The holes are assumed to be unresolved. They therefore cover one pixel and have
        value zero, i.e. no background emission.
        """
        hdu = fits.ImageHDU()
        hdu.header.update(header)
        hdu.data = np.ones((2047, 2047)) * self.pixarea.value # TODO test - do we need to go back to 2048?
        tab = self.data_container.table
        xpix = (tab['x'] - header['CRVAL1']) / header['CDELT1'] + header['CRPIX1'] - 1
        ypix = (tab['y'] - header['CRVAL2']) / header['CDELT2'] + header['CRPIX2'] - 1
        hdu.data[ypix.astype(int), xpix.astype(int)] = 0
        return hdu
