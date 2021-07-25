import warnings
from copy import deepcopy

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

from scopesim import utils, rc
from scopesim.optics import image_plane_utils as imp_utils


def is_field_in_fov(fov_header, table_or_imagehdu, wcs_suffix=""):
    """
    Returns True if Source.field footprint is inside the FieldOfView footprint

    Parameters
    ----------
    fov_header : fits.Header
        Header from a FieldOfView object
    table_or_imagehdu : [astropy.Table, astropy.ImageHDU]
        Field object from a Source object
    wcs_suffix : str
        ["S", "D"] Coordinate system: Sky or Detector

    Returns
    -------
    is_inside_fov : bool

    """

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

    image = np.zeros((fov_header["NAXIS2"], fov_header["NAXIS1"]))
    canvas_hdu = fits.ImageHDU(header=fov_header, data=image)
    order = int(rc.__config__["!SIM.computing.spline_order"])
    pixel_area = fov_header["CDELT1"] * fov_header["CDELT2"] * \
                 u.Unit(fov_header["CUNIT1"]).to(u.arcsec) ** 2

    for ii in fields_indexes:
        field = src.fields[ii]
        if isinstance(field, fits.ImageHDU):
            ref = field.header["SPEC_REF"]
            flux = src.photons_in_range(wave_min, wave_max, area, indexes=[ref])
            image = np.zeros((fov_header["NAXIS2"], fov_header["NAXIS1"]))
            temp_hdu = fits.ImageHDU(header=fov_header, data=image)

            if field.header.get("BG_SRC", False) and \
                    field.header["NAXIS1"] <= 1 and \
                    field.header["NAXIS2"] <= 1:
                # .. todo: check if we need to take pixel_scale into account
                temp_hdu.data += flux[0].value * pixel_area
            else:
                temp_hdu = imp_utils.add_imagehdu_to_imagehdu(field, temp_hdu,
                                                              order, wcs_suffix)
                temp_hdu.data *= flux[0].value

            canvas_hdu.data += temp_hdu.data

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


def extract_common_field(field, fov_volume):
    """
    Extracts the overlapping parts of a field within a FOV volume

    Parameters
    ----------
    field : Table or ImageHDU
    fov_volume : dict
        Contains {"xs": [xmin, xmax], "ys": [ymin, ymax],
                  "waves": [wave_min, wave_max],
                  "xy_unit": "deg" or "mm", "wave_unit": "um"}

    Returns
    -------
    field_new : Table or ImageHDU

    """
    if isinstance(field, Table):
        mask = (field["x"] >= fov_volume["xs"][0]) * \
               (field["x"] < fov_volume["xs"][1]) * \
               (field["y"] >= fov_volume["ys"][0]) * \
               (field["y"] < fov_volume["ys"][1])
        field_new = field[mask]
    elif isinstance(field, fits.ImageHDU):
        field_new = extract_area_from_imagehdu(field, fov_volume)
    else:
        raise ValueError("field must be either Table or ImageHDU: {}"
                         "".format(type(field)))

    return field_new


def extract_area_from_imagehdu(imagehdu, fov_volume):
    """
    Extracts the part of a ImageHDU that fits inside the fov_volume

    Parameters
    ----------
    imagehdu : fits.ImageHDU
        The field ImageHDU, either an image of a wavelength [um] cube
    fov_volume : dict
        Contains {"xs": [xmin, xmax], "ys": [ymin, ymax],
                  "waves": [wave_min, wave_max],
                  "xy_unit": "deg" or "mm", "wave_unit": "um"}

    Returns
    -------
    new_imagehdu : fits.ImageHDU

    """
    hdr = imagehdu.header
    new_hdr = {}

    x_hdu, y_hdu = imp_utils.calc_footprint(imagehdu)  # field edges in "deg"
    x_fov, y_fov = fov_volume["xs"], fov_volume["ys"]

    x0s, x1s = max(min(x_hdu), min(x_fov)), min(max(x_hdu), max(x_fov))
    y0s, y1s = max(min(y_hdu), min(y_fov)), min(max(y_hdu), max(y_fov))

    xp, yp = imp_utils.val2pix(hdr, np.array([x0s, x1s]), np.array([y0s, y1s]))
    (x0p, x1p), (y0p, y1p) = np.round(xp).astype(int), np.round(yp).astype(int)

    new_hdr = imp_utils.header_from_list_of_xy([x0s, x1s], [y0s, y1s],
                                               pixel_scale=hdr["CDELT1"])

    if hdr["NAXIS"] == 3:
        wval, wdel, wpix, wlen = [hdr[kw] for kw in ["CRVAL3", "CDELT3",
                                                     "CRPIX3", "NAXIS3"]]
        # ASSUMPTION - cube wavelength is in regularly spaced units of um
        w_hdu = [wval - wdel * wpix, wval + wdel * (wlen - wpix)]
        w_fov = fov_volume["waves"]

        w0s, w1s = max(min(w_hdu), min(w_fov)), min(max(w_hdu), max(w_fov))
        wp = (np.array([w0s, w1s]) - wval) / wdel + wpix
        (w0p, w1p) = np.round(wp).astype(int)

        new_hdr.update({"NAXIS": 3,
                        "NAXIS3": w1p - w0p,
                        "CRVAL3": w0s,
                        "CRPIX3": 0,
                        "CDELT3": hdr["CDELT3"],
                        "CUNIT3": hdr["CUNIT3"],
                        "CTYPE3": hdr["CTYPE3"]})

        data = imagehdu.data[w0p:w1p, y0p:y1p, x0p:x1p]
    else:
        data = imagehdu.data[y0p:y1p, x0p:x1p]

    new_imagehdu = fits.ImageHDU(data=data)
    new_imagehdu.header.update(new_hdr)

    return new_imagehdu

