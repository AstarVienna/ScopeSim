import numpy as np


from astropy.io import fits
from astropy.table import Table

from .image_plane_utils import add_table_to_imagehdu, add_imagehdu_to_imagehdu

from ..base_classes import ImagePlaneBase
from .. import rc
from .. import utils


class ImagePlane(ImagePlaneBase):
    """
    A class to act as a canvas onto which to project `Source` images or tables

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

    def __init__(self, header, **kwargs):

        max_seg_size = rc.__config__["!SIM.computing.max_segment_size"]
        self.meta = {"SIM_MAX_SEGMENT_SIZE" : max_seg_size}
        self.meta.update(kwargs)
        self.id = header["IMGPLANE"] if "IMGPLANE" in header else 0

        if not any([utils.has_needed_keywords(header, s)
                    for s in ["", "D", "S"]]):
            raise ValueError("header must have a valid image-plane WCS: {}"
                             "".format(dict(header)))

        image = np.zeros((header["NAXIS2"]+1, header["NAXIS1"]+1))
        self.hdu = fits.ImageHDU(data=image, header=header)

    def add(self, hdus_or_tables, sub_pixel=None, spline_order=None, wcs_suffix=""):
        """
        Add a projection of an image or table files to the canvas

        .. note::
            If a Table is provided, it must include the following columns:
            `x_mm`, `y_mm`, and `flux`.

            Units for the columns should be provided in the
            <Table>.unit attribute or as an entry in the table's meta dictionary
            using this syntax: <Table>.meta["<colname>_unit"] = <unit>.

            For example::

              tbl["x"].unit = u.arcsec   # or
              tbl.meta[x_unit"] = "deg"

            If no units are given, default units will be assumed. These are:

            - `x`, `y`: `arcsec`
            - `flux` : `ph / s / pix`

        Parameters
        ----------
        hdus_or_tables : `fits.ImageHDU` or `astropy.Table`
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
            sub_pixel = utils.from_currsys("!SIM.sub_pixel.flag")
        if spline_order is None:
            spline_order = utils.from_currsys("!SIM.computing.spline_order")

        if isinstance(hdus_or_tables, (list, tuple)):
            for hdu_or_table in hdus_or_tables:
                self.add(hdu_or_table, sub_pixel, spline_order, wcs_suffix)
        else:
            if isinstance(hdus_or_tables, Table):
                self.hdu.header["COMMENT"] = "Adding files from table"
                self.hdu = add_table_to_imagehdu(hdus_or_tables, self.hdu,
                                                 sub_pixel, wcs_suffix)
            elif isinstance(hdus_or_tables, fits.ImageHDU):
                self.hdu.header["COMMENT"] = "Adding files from table"
                self.hdu = add_imagehdu_to_imagehdu(hdus_or_tables, self.hdu,
                                                    spline_order, wcs_suffix)

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
