# -*- coding: utf-8 -*-
"""Contains ``ImagePlane`` class."""

from warnings import warn

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from .image_plane_utils import add_table_to_imagehdu, add_imagehdu_to_imagehdu

from ..utils import (
    from_currsys,
    has_needed_keywords,
    get_logger,
    zeros_from_header,
    image_plotter,
    cube_plotter,
)


logger = get_logger(__name__)


class ImagePlane:
    """
    A class to act as a canvas onto which to project `Source` images or tables.

    Parameters
    ----------
    header : `fits.Header`
        Must contain a valid WCS

    .. todo: Write the code to deal with a canvas larger than max_segment_size

    Examples
    --------
    ::

        from astropy.table import Table
        from scopesim.optics import image_plane as imp

        my_point_source_table = Table(names=["x", "y", "flux"],
                                      data=[(0,  1,  2)*u.mm,
                                            (0, -5, 10)*u.mm,
                                            (100,50,25)*u.ph/u.s])

        hdr = imp.make_image_plane_header([my_point_source_table],
                                           pixel_size=0.015*u.mm)
        img_plane = imp.ImagePlane(hdr)
        img_plane.add(my_point_source_table)

        print(img_plane.image)

    """

    def __init__(self, header, cmds=None, **kwargs):

        self.cmds = cmds
        self.meta = {}
        self.meta.update(kwargs)
        self.id = header.get("IMGPLANE", 0)

        if not any(has_needed_keywords(header, s)
                   for s in ["", "D", "S"]):
            raise ValueError(f"header must have a valid image-plane WCS: "
                             f"{dict(header)}")

        image = zeros_from_header(header)
        self.hdu = fits.ImageHDU(data=image, header=header)
        self.hdu.header["BUNIT"] = "ph s-1"  # photons per second (per pixel)

        self._det_wcs = self._get_wcs(header, "D")
        logger.debug("det %s", self._det_wcs)
        self._sky_wcs = self._get_wcs(header, " ")
        logger.debug("sky %s", self._sky_wcs)

    def add(self, hdus, sub_pixel=None, spline_order=None,
            wcs_suffix=""):
        """
        Add a projection of an image to the canvas.

        .. versionchanged:: 0.10.0

           Adding a table directly to the ImagePlane is deprecated. Use FOV to
           add tables and image HDUs together before adding them to here.

        .. versionchanged:: 0.12.0

           Support for adding tables fully removed. See note above.

        Parameters
        ----------
        hdus : `fits.ImageHDU` or list thereof
            The input to be projected onto the image plane. See above.

        sub_pixel : bool, optional
            Default is False. Dictates if point files should be projected with
            sub-pixel shifts or not. Accounting for sub-pixel shifts is approx.
            5x slower.

        spline_order : int, optional
            Order of spline interpolations used in ``scipy.ndimage`` functions
            ``zoom`` and ``rotate``.

        wcs_suffix : str, optional
            Default "". For sky coords - "" or "S", Detector coords - "D"

        """
        if sub_pixel is None:
            sub_pixel = from_currsys("!SIM.sub_pixel.flag", self.cmds)
        if spline_order is None:
            spline_order = from_currsys("!SIM.computing.spline_order", self.cmds)

        if isinstance(hdus, (list, tuple)):
            logger.debug("Adding multiple HDUs to ImagePlane.")
            for hdu in hdus:
                self.add(hdu, sub_pixel, spline_order, wcs_suffix)
        else:
            if not isinstance(hdus, fits.ImageHDU):
                raise TypeError("Only ImageHDUs may be added to ImagePlane.")

            logger.debug("Adding HDU with shape %d to ImagePlane.",
                         hdus.data.shape)
            self.hdu = add_imagehdu_to_imagehdu(
                hdus, self.hdu, spline_order, wcs_suffix)
            if (img_bunit := hdus.header.get("BUNIT")) is not None:
                if img_bunit != (imp_bunit := self.hdu.header["BUNIT"]):
                    logger.warning("Added mismatched BUNIT %s to %s.",
                                   img_bunit, imp_bunit)
            else:
                logger.info("No BUNIT found in added HDU.")

    @property
    def header(self):
        return self.hdu.header

    @property
    def data(self):
        return self.hdu.data

    @data.setter
    def data(self, new_data):
        self.hdu.data = new_data

    @property
    def image(self):
        return self.data

    def view(self, sub_pixel):
        # for consistency with FieldOfView
        return self.data

    @staticmethod
    def _get_wcs(header: fits.Header, key: str) -> WCS:
        sky_alias = {" ", "S"}
        try:
            wcs = WCS(header, key=key)
        except KeyError:
            # retry with alias
            sky_alias.discard(key)
            try:
                wcs = WCS(header, key=sky_alias.pop())
            except KeyError:
                wcs = None
        return wcs

    def plot(self):
        """Plot data in image plane.

        .. versionadded:: 0.11.0

        """
        if self.header["NAXIS"] == 3:
            return cube_plotter(self.hdu)
        return image_plotter(self.hdu)
