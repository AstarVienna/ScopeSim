import warnings

import numpy as np

from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

from . import image_plane_utils as imp_utils
from ..source.source2 import Source

from .. import utils
from .. import rc


class FieldOfView:
    """
    A FOV is a monochromatic image. Flux units after extracting the fields from
    the Source are in ph/s/pixel

    The initial header should contain an on-sky WCS description:
    - CDELT-, CUNIT-, NAXIS- : for pixel scale and size (assumed CUNIT in deg)
    - CRVAL-, CRPIX- : for positioning the final image
    - CTYPE- : is assumed to be "RA---TAN", "DEC---TAN"

    and an image-plane WCS description
    - CDELT-D, CUNIT-D, NAXISn : for pixel scale and size (assumed CUNIT in mm)
    - CRVAL-D, CRPIX-D : for positioning the final image
    - CTYPE-D : is assumed to be "LINEAR", "LINEAR"

    The wavelength range is given by waverange

    """

    def __init__(self, header, waverange, **kwargs):
        self.meta = {"id" : None,
                     "wave_binwidth" : rc.__rc__["SIM_SPEC_RESOLUTION"],
                     "wave_min" : utils.quantify(waverange[0], u.um),
                     "wave_max" : utils.quantify(waverange[1], u.um),
                     "area" : 1 * u.m**2,
                     "sub_pixel" : rc.__rc__["SIM_SUB_PIXEL_ACCURACY"]}
        self.meta.update(kwargs)

        if not any([utils.has_needed_keywords(header, s) for s in ["", "S"]]):
            raise ValueError("header must contain a valid sky-plane WCS: {}"
                             "".format(dict(header)))
        if not utils.has_needed_keywords(header, "D"):
            raise ValueError("header must contain a valid image-plane WCS: {}"
                             "".format(dict(header)))

        self.hdu = fits.ImageHDU(header=header)
        self.hdu.header["NAXIS"] = 2
        self.hdu.header["NAXIS1"] = header["NAXIS1"]
        self.hdu.header["NAXIS2"] = header["NAXIS2"]

        self.fields = []

    def extract_from(self, src):
        """ ..assumption: Bandpass has been applied"""
        
        if not isinstance(src, Source):
            raise ValueError("source must be a Source object: {}"
                             "".format(type(src)))

        wave_min = utils.quantify(self.meta["wave_min"], u.um).value
        wave_max = utils.quantify(self.meta["wave_max"], u.um).value
        area = self.meta["area"]

        # determine which fields are inside the field of view
        fields_mask = [is_field_in_fov(self.hdu.header, field)
                       for field in src.fields]
        fields_indexes = np.where(fields_mask)[0]
        tbl_fields_mask = np.array([isinstance(field, Table)
                                    for field in src.fields])
        img_fields_mask = np.array([isinstance(field, fits.ImageHDU)
                                    for field in src.fields])

        # combine all Table fields
        if sum(tbl_fields_mask * fields_mask) > 0:
            combined_table = combine_table_fields(self.hdu.header, src,
                                                  fields_indexes)
            tbl = make_flux_table(combined_table, src, wave_min, wave_max, area)
            xd, yd = sky2fp(self.hdu.header, tbl["x"], tbl["y"])
            tbl.add_columns([Column(name="x_mm", data=xd, unit=u.mm),
                             Column(name="y_mm", data=yd, unit=u.mm)])
            self.fields += [tbl]

        # combine all ImageHDU fields
        if sum(img_fields_mask * fields_mask) > 0:
            imagehdu = combine_imagehdu_fields(self.hdu.header, src,
                                               fields_indexes, wave_min,
                                               wave_max, area)
            self.fields += [imagehdu]

    def view(self, sub_pixel=False):
        self.hdu.data = np.zeros((self.hdu.header["NAXIS1"],
                                  self.hdu.header["NAXIS2"]))
        if len(self.fields) > 0:
            for field in self.fields:
                if isinstance(field, Table):
                    self.hdu = imp_utils.add_table_to_imagehdu(field, self.hdu,
                                                               sub_pixel)
                elif isinstance(field, fits.ImageHDU):
                    self.hdu.data += field.data

        return self.hdu.data

    @property
    def header(self):
        return self.hdu.header

    @property
    def data(self):
        if self.hdu.data is None:
            self.view(self.meta["sub_pixel"])
        return self.hdu.data

    @property
    def image(self):
        return self.data

    @property
    def corners(self):
        sky_corners = imp_utils.calc_footprint(self.header)
        imp_corners = imp_utils.calc_footprint(self.header, "D")
        return sky_corners, imp_corners

    def __repr__(self):
        msg = "FOV id: {}, with dimensions ({}, {})\n" \
              "".format(self.meta["id"], self.header["NAXIS1"],
                        self.header["NAXIS2"])
        msg += "Sky centre: ({},{})\n" \
               "".format(self.header["CRVAL1"], self.header["CRVAL2"])
        msg += "Image centre: ({},{})\n" \
               "".format(self.header["CRVAL1D"], self.header["CRVAL2D"])

        return msg


def is_field_in_fov(fov_header, table_or_imagehdu, wcs_suffix=""):

    s = wcs_suffix
    pixel_scale = utils.quantify(fov_header["CDELT1"+s], u.deg)

    if isinstance(table_or_imagehdu, Table):
        ext_hdr = imp_utils._make_bounding_header_for_tables(
                                            [table_or_imagehdu], pixel_scale)
    elif isinstance(table_or_imagehdu, fits.ImageHDU):
        ext_hdr = imp_utils._make_bounding_header_from_imagehdus(
                                            [table_or_imagehdu], pixel_scale)
    else:
        warnings.warn("Input was neither Table nor ImageHDU: {}"
                      "".format(table_or_imagehdu))
        return False

    ext_xsky, ext_ysky = imp_utils.calc_footprint(ext_hdr, wcs_suffix)
    fov_xsky, fov_ysky = imp_utils.calc_footprint(fov_header, wcs_suffix)

    is_inside_fov = min(ext_xsky) < max(fov_xsky) and \
                    max(ext_xsky) > min(fov_xsky) and \
                    min(ext_ysky) < max(fov_ysky) and \
                    max(ext_ysky) > min(fov_ysky)

    return is_inside_fov


def make_flux_table(source_tbl, src, wave_min, wave_max, area):
    fluxes = np.zeros(len(src.spectra))
    ref_set = list(set(source_tbl["ref"]))
    flux_set = src.photons_in_range(wave_min, wave_max, area, ref_set)
    fluxes[ref_set] = flux_set

    ref = source_tbl["ref"]
    weight = source_tbl["weight"]
    flux_col = Column(name="flux", data=fluxes[ref] * weight)
    x_col = source_tbl["x"]
    y_col = source_tbl["y"]

    tbl = Table()
    tbl.add_columns([x_col, y_col, flux_col])

    return tbl


def combine_table_fields(fov_header, src, field_indexes):
    """
    Combines a list of Table objects into a single one bounded by the Header WCS

    Parameters
    ----------
    fov_header : fits.Header
        Header from a FieldOfView objects
    src : Source object
    field_indexes : list of int

    Returns
    -------
    tbl : Table

    """

    fov_xsky, fov_ysky = imp_utils.calc_footprint(fov_header)

    x, y, ref, weight = [], [], [], []

    for ii in field_indexes:
        field = src.fields[ii]
        if isinstance(field, Table):
            xcol = utils.quantity_from_table("x", field, u.arcsec)
            ycol = utils.quantity_from_table("y", field, u.arcsec)
            x += list(xcol.to(u.deg).value)
            y += list(ycol.to(u.deg).value)
            ref += list(field["ref"])
            weight += list(field["weight"])

    x = np.array(x)
    y = np.array(y)
    mask = np.array(x < max(fov_xsky)) * np.array(x > min(fov_xsky)) * \
           np.array(y < max(fov_ysky)) * np.array(y > min(fov_ysky))

    x = x[mask]
    y = y[mask]
    ref = np.array(ref)[mask]
    weight = np.array(weight)[mask]

    tbl = Table(names=["x", "y", "ref", "weight"], data=[x, y, ref, weight])
    tbl["x"].unit = u.deg
    tbl["y"].unit = u.deg

    return tbl


def combine_imagehdu_fields(fov_header, src, fields_indexes, wave_min, wave_max,
                            area, wcs_suffix=""):
    """
    Combines a list of ImageHDUs into a single one bounded by the Header WCS

    Parameters
    ----------
    fov_header : fits.Header
        Header from the FieldOfView
    src : Source object
    fields_indexes : list of ints
        Which indexes from <Source>.fields to use
    wave_min : float
        [deg] Blue spectral border
    wave_max : float
        [deg] Red spectral border
    area : float
        [m2] Area of the primary aperture
    wcs_suffix : str
        Which coordinate system to use
        - "" for the on-sky coordinate system
        - "D" for the image-plane coordinate system

    Returns
    -------
    canvas_hdu : fits.ImageHDU

    """

    image = np.zeros((fov_header["NAXIS1"], fov_header["NAXIS2"]))
    canvas_hdu = fits.ImageHDU(header=fov_header, data=image)
    order = int(rc.__rc__["SIM_SPLINE_ORDER"])

    for ii in fields_indexes:
        if isinstance(src.fields[ii], fits.ImageHDU):
            ref = src.fields[ii].header["SPEC_REF"]
            flux = src.photons_in_range(wave_min, wave_max, area, indexes=[ref])
            image = np.zeros((fov_header["NAXIS1"], fov_header["NAXIS2"]))
            temp_hdu = fits.ImageHDU(header=fov_header, data=image)
            temp_hdu = imp_utils.add_imagehdu_to_imagehdu(src.fields[ii],
                                                          temp_hdu, order=order,
                                                          wcs_suffix=wcs_suffix)
            canvas_hdu.data += temp_hdu.data * flux[0].value

    return canvas_hdu


def sky2fp(header, xsky, ysky):
    """
    Convert sky coordinates to image plane coordinated

    Parameters
    ----------
    header : Header
        Header of a FieldOfView object which contains two sets of WCS keywords
    xsky, ysky : float, array
        [deg] The on-sky coordinated

    Returns
    -------
    xdet, ydet : float, array
        [mm] The coordinated on the image plane

    """

    xpix, ypix = imp_utils.val2pix(header, xsky, ysky)
    xdet, ydet = imp_utils.pix2val(header, xpix, ypix, "D")

    return xdet, ydet
