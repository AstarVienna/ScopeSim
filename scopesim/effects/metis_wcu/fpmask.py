"""A class for the METIS WCU focal-plane mask"""

from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from ..data_container import DataContainer
from ...utils import find_file, from_currsys, get_logger, figure_factory
from ...optics.image_plane_utils import sub_pixel_fractions

logger = get_logger(__name__)

class FPMask:
    """Focal-plane mask for the METIS WCU

    Parameters
    ----------
    See :class:`DataContainer` for input parameters

    """
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

    def __init__(self,
                 maskname: Path | str | None = None,
                 fpmask_filename_format: str | None = None,
                 angle: float = 0,
                 shift: tuple = (0, 0),
                 **kwargs
                 ):
        logger.debug("Initialising FPMask with %s", maskname)
        self.name = maskname
        self.angle = angle
        self.shift = shift
        self.xpix = []
        self.ypix = []
        if maskname == "open":
            self.holehdu = fits.ImageHDU()
            self.holehdu.header.update(self.hdr)
            self.opaquehdu = None
        else:
            # Try to find the file as a path
            if find_file(maskname, silent=True) is None:
                file_format = from_currsys(fpmask_filename_format)
                self.filename = file_format.format(maskname)
            else:
                self.filename = maskname

            self.data_container = DataContainer(filename=self.filename, **kwargs)
            self.pixarea = (self.hdr['CDELT1'] * u.Unit(self.hdr['CUNIT1'])
                            * self.hdr['CDELT2'] * u.Unit(self.hdr['CUNIT2']))
            self.make_hdus(header=self.hdr)


    def make_hdus(self, header):
        """Create an hdu for the holes in fpmask

        The holes are assumed to be unresolved. They therefore cover one pixel and have
        a value corresponding to the actual solid angle covered by the hole.
        """
        holehdu = fits.ImageHDU()
        holehdu.header.update(header)
        holehdu.data = np.zeros((2047, 2047))

        opaquehdu = fits.ImageHDU()
        opaquehdu.header.update(header)
        opaquehdu.data = np.ones((2047, 2047)) * self.pixarea.value

        # Hole locations
        tab = self.data_container.table
        xhole = tab['x'].data
        yhole = tab['y'].data
        diam = tab['diam'].data

        if self.angle != 0:
            rangle = np.deg2rad(self.angle)
            xtmp = xhole * np.cos(rangle) - yhole * np.sin(rangle)
            ytmp = xhole * np.sin(rangle) + yhole * np.cos(rangle)
            xhole = xtmp
            yhole = ytmp

        xhole += self.shift[0]
        yhole += self.shift[1]

        xpix = (xhole - header['CRVAL1']) / header['CDELT1'] + header['CRPIX1'] - 1
        ypix = (yhole - header['CRVAL2']) / header['CDELT2'] + header['CRPIX2'] - 1
        in_field = (xpix > 0) * (xpix < 2047) * (ypix > 0) * (ypix < 2047)


        for x, y, d in zip(xpix[in_field], ypix[in_field], diam[in_field]):
            holearea = (d/2)**2 * np.pi
            xint, yint, fracs = sub_pixel_fractions(x, y)
            holehdu.data[yint, xint] = np.array(fracs) * holearea
            opaquehdu.data[yint, xint] = 0
        self.xpix = xpix
        self.ypix = ypix
        self.holehdu = holehdu
        self.opaquehdu = opaquehdu


    def plot(self):
        """Plot the location of the holes"""
        _, axes = figure_factory()
        axes.plot(self.xpix, self.ypix, 'o')
        axes.set_xlim(0, 2048)
        axes.set_ylim(0, 2048)
        axes.set_aspect('equal')

        return axes

    def __str__(self) -> str:
        return f"""{self.__class__.__name__}: "{self.name}"
    Angle:        {self.angle} deg
    Shift:        {self.shift} arcsec"""
