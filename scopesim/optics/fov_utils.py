# -*- coding: utf-8 -*-

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

from . import image_plane_utils as imp_utils
from ..utils import from_currsys, quantity_from_table


# TODO: Unused function. Remove?
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


# TODO: Unused function. Remove?
def combine_table_fields(fov_header, src, field_indexes):
    """
    Combine list of Table objects into a single one bounded by the Header WCS.

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
    xy = imp_utils.calc_footprint(fov_header)
    fov_xsky, fov_ysky = xy[:, 0], xy[:, 1]

    x, y, ref, weight = [], [], [], []

    for ii in field_indexes:
        field = src.fields[ii]
        if isinstance(field, Table):
            xcol = quantity_from_table("x", field, u.arcsec)
            ycol = quantity_from_table("y", field, u.arcsec)
            x += list(xcol.to(u.deg).value)
            y += list(ycol.to(u.deg).value)
            ref += list(field["ref"])
            weight += list(field["weight"])

    x = np.array(x)
    y = np.array(y)
    mask = (np.array(x < max(fov_xsky)) * np.array(x > min(fov_xsky)) *
            np.array(y < max(fov_ysky)) * np.array(y > min(fov_ysky)))

    x = x[mask]
    y = y[mask]
    ref = np.array(ref)[mask]
    weight = np.array(weight)[mask]

    tbl = Table(names=["x", "y", "ref", "weight"], data=[x, y, ref, weight])
    tbl["x"].unit = u.deg
    tbl["y"].unit = u.deg

    return tbl


# TODO: Unused functionv (except in a comment somewhere). Remove?
def combine_imagehdu_fields(fov_header, src, fields_indexes, wave_min, wave_max,
                            area, wcs_suffix="", cmds=None):
    """
    Combine list of ImageHDUs into a single one bounded by the Header WCS.

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
    spline_order = from_currsys("!SIM.computing.spline_order", cmds)
    pixel_area = (fov_header["CDELT1"] * fov_header["CDELT2"] *
                  u.Unit(fov_header["CUNIT1"].lower()).to(u.arcsec) ** 2)

    for ii in fields_indexes:
        field = src.fields[ii]
        if isinstance(field, fits.ImageHDU):
            ref = field.header["SPEC_REF"]
            flux = src.photons_in_range(wave_min, wave_max, area,
                                        indexes=[ref])
            image = np.zeros((fov_header["NAXIS2"], fov_header["NAXIS1"]))
            temp_hdu = fits.ImageHDU(header=fov_header, data=image)

            if field.header.get("BG_SRC", False) and \
                    field.header["NAXIS1"] <= 1 and \
                    field.header["NAXIS2"] <= 1:
                # .. todo: check if we need to take pixel_scale into account
                temp_hdu.data += flux[0].value * pixel_area
            else:
                temp_hdu = imp_utils.add_imagehdu_to_imagehdu(
                    field, temp_hdu, spline_order, wcs_suffix)
                temp_hdu.data *= flux[0].value

            canvas_hdu.data += temp_hdu.data

    return canvas_hdu
