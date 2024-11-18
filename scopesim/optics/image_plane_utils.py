# -*- coding: utf-8 -*-

from itertools import product
from collections.abc import Iterable

import numpy as np
from astropy import units as u
from astropy.wcs import WCS, find_all_wcs
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from scipy import ndimage as ndi

from ..utils import (unit_from_table, quantity_from_table, has_needed_keywords,
                     get_logger)


logger = get_logger(__name__)


def get_canvas_header(hdu_or_table_list, pixel_scale=1 * u.arcsec):
    """
    Generate a fits.Header with a WCS that covers everything in the FOV.

    Parameters
    ----------
    hdu_or_table_list : list
        A list of Tables and/or ImageHDU py_objects

    pixel_scale : astropy.Quantity
        [arcsec] The pixel scale of the projection. Default in 1 arcsec

    Returns
    -------
    header : fits.Header
        A Header containing a WCS and NAXISn values to build an ImageHDU

    """
    size_warning = ("Header dimension are {adverb} large: {num_pix}. "
                    "Any image made from this header will use more that "
                    ">{size} in memory")

    def _get_headers(hdus_or_tables):
        tables = []
        for hdu_or_table in hdus_or_tables:
            if isinstance(hdu_or_table, fits.ImageHDU):
                yield hdu_or_table.header
            elif isinstance(hdu_or_table, fits.Header):
                yield hdu_or_table
            elif isinstance(hdu_or_table, Table):
                tables.append(hdu_or_table)
            else:
                raise TypeError(
                    "hdu_or_table_list may only contain fits.ImageHDU, Table "
                    f"or fits.Header, found {type(hdu_or_table)}.")
        if tables:
            yield _make_bounding_header_for_tables(*tables,
                                                   pixel_scale=pixel_scale)

    headers = list(_get_headers(hdu_or_table_list))

    if not headers:
        logger.warning("No tables or ImageHDUs were passed")
        return None

    hdr = _make_bounding_header_from_headers(*headers, pixel_scale=pixel_scale)

    num_pix = hdr["NAXIS1"] * hdr["NAXIS2"]
    if num_pix > 2 ** 28:
        raise MemoryError(size_warning.format(adverb="too", num_pix=num_pix,
                                              size="8 GB"))
    if num_pix > 2 ** 25:  # 2 * 4096**2
        logger.warning(size_warning.format(adverb="", num_pix=num_pix,
                                           size="256 MB"))
    return hdr


def _make_bounding_header_from_headers(*headers, pixel_scale=1*u.arcsec):
    """
    Return a Header with WCS and NAXISn keywords bounding all input ImageHDUs.

    Parameters
    ----------
    headers : list of fits.ImageHDU
    pixel_scale : u.Quantity
        [arcsec]

    Returns
    -------
    hdr : fits.Header

    """
    wcs_suffix = "D" if pixel_scale.unit.physical_type == "length" else ""
    unit = u.Unit(_get_unit_from_headers(*headers, wcs_suffix=wcs_suffix))

    if unit.physical_type == "angle":
        unit = "deg"
        pixel_scale = pixel_scale.to(u.deg).value
    else:
        pixel_scale = pixel_scale.to(unit).value

    extents = [calc_footprint(header, wcs_suffix, unit) for header in headers]
    pnts = np.vstack(extents)

    hdr = header_from_list_of_xy(pnts[:, 0], pnts[:, 1],
                                 pixel_scale, wcs_suffix)
    hdr["NAXIS1"] += 1
    hdr["NAXIS2"] += 1
    hdr[f"CRVAL1{wcs_suffix}"] -= 0.5 * pixel_scale
    hdr[f"CRVAL2{wcs_suffix}"] -= 0.5 * pixel_scale

    return hdr


def _make_bounding_header_for_tables(*tables, pixel_scale=1*u.arcsec):
    """
    Return a Header with WCS and NAXISn keywords bounding all input Tables.

    Parameters
    ----------
    tables : list of astropy.Tables
        [arcsec] Must contain columns: "x", "y"

    pixel_scale : u.Quantity
        [arcsec]

    Returns
    -------
    hdr : fits.Header

    """
    wcs_suffix = "D" if pixel_scale.unit.physical_type == "length" else ""
    # FIXME: Convert to deg here? If yes, remove the arcsec=True below...
    # Note: this could all be a lot simpler if we have consistent units, i.e.
    #       don't need to convert mm -> mm and arcsec -> arcsec (or deg)
    new_unit = u.mm if wcs_suffix == "D" else u.arcsec  # u.deg
    tbl_unit = u.mm if wcs_suffix == "D" else u.arcsec
    x_name = "x_mm" if wcs_suffix == "D" else "x"
    y_name = "y_mm" if wcs_suffix == "D" else "y"

    pixel_scale = pixel_scale.to(new_unit)

    extents = []
    for table in tables:
        extent = calc_table_footprint(table, x_name, y_name,
                                      tbl_unit, new_unit,
                                      padding=pixel_scale)
        extents.append(extent)
    pnts = np.vstack(extents)

    # TODO: check if this could just use create_wcs_from_points
    hdr = header_from_list_of_xy(pnts[:, 0], pnts[:, 1], pixel_scale.value,
                                 wcs_suffix, arcsec=wcs_suffix != "D")
    return hdr


def create_wcs_from_points(points: np.ndarray,
                           pixel_scale: float,
                           wcs_suffix: str = "") -> tuple[WCS, np.ndarray]:
    """
    Create `astropy.wcs.WCS` instance that fits all points inside.

    Parameters
    ----------
    corners : (N, 2) array
        2D array of N >= 2 points in the form of [x, y].
    pixel_scale : float
        DESCRIPTION.
    wcs_suffix : str, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    new_wcs : TYPE
        Newly created WCS instance.
    naxis : TYPE
        Array of NAXIS needed to fit all points.

    """
    # TODO: should pixel_scale be called pixel_size by conventions elsewhere?
    # TODO: add quantity stuff to docstring
    # TODO: find out how this plays with chunks
    if wcs_suffix != "D":
        points = _fix_360(points)

    # TODO: test whether abs(pixel_scale) breaks anything
    pixel_scale = abs(pixel_scale)
    extent = points.ptp(axis=0) / pixel_scale
    naxis = extent.round().astype(int)

    # FIXME: Woule be nice to have D headers referenced at (1, 1), but that
    #        breaks some things, likely down to rescaling...
    # if wcs_suffix == "D":
    #     crpix = np.array([1., 1.])
    #     crval = points.min(axis=0)
    # else:

    if isinstance(naxis, u.Quantity):
        naxis = naxis.decompose()
        if naxis.unit != "pixel":
            raise u.UnitConversionError("If given quantities, must resolve to "
                                        f"pix, got '{naxis.unit}' instead.")
        naxis = naxis.value

    crpix = (naxis + 1) / 2
    crval = (points.min(axis=0) + points.max(axis=0)) / 2

    ctype = "LINEAR" if wcs_suffix in "DX" else "RA---TAN"

    if isinstance(points, u.Quantity):
        cunit = points.unit
    else:
        if wcs_suffix == "D":
            cunit = "mm"
        else:
            cunit = "deg"

    if isinstance(pixel_scale, u.Quantity):
        pixel_scale = pixel_scale.value

    new_wcs = WCS(key=wcs_suffix)
    new_wcs.wcs.ctype = 2 * [ctype]
    new_wcs.wcs.cunit = 2 * [cunit]
    new_wcs.wcs.cdelt = np.array(2 * [pixel_scale])
    new_wcs.wcs.crval = crval
    new_wcs.wcs.crpix = crpix

    return new_wcs, naxis


def header_from_list_of_xy(x, y, pixel_scale, wcs_suffix="", arcsec=False):
    """
    Make a header large enough to contain all x,y on-sky coordinates.

    Parameters
    ----------
    x, y : list of floats
        [deg, mm] List of sky coordinates to be bounded by the NAXISn keys
    pixel_scale : float
        [deg, mm]

    Returns
    -------
    hdr : fits.Header

    """
    points = np.column_stack((x, y))

    if arcsec:
        points <<= u.arcsec
        pixel_scale <<= u.arcsec / u.pixel

    new_wcs, naxis = create_wcs_from_points(points, pixel_scale, wcs_suffix)

    hdr = fits.Header()
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = int(naxis[0])
    hdr["NAXIS2"] = int(naxis[1])
    hdr.update(new_wcs.to_header())

    return hdr


def add_table_to_imagehdu(table: Table, canvas_hdu: fits.ImageHDU,
                          sub_pixel: bool = True,
                          wcs_suffix: str = "") -> fits.ImageHDU:
    """
    Add files from an astropy.Table to the image of an fits.ImageHDU.

    Parameters
    ----------
    table : astropy.Table
        Must contain the columns "x_mm", "y_mm", "flux" with the units in the
        column attribute .unit, or in the table.meta dictionary as
        "<colname>_unit". Default units are ``mm`` and ``ph / s / pix``

    canvas_hdu : fits.ImageHDU
        The ImageHDU onto which the table files should be projected.
        This must include a valid WCS

    sub_pixel : bool, optional
        Default is True. If True, sub-pixel shifts of files will be taken into
        account when projecting onto the canvas pixel grid. This takes about 5x
        longer than ignoring the sub-pixel shifts

    wcs_suffix : str, optional

    Returns
    -------
    canvas_hdu : fits.ImageHDU

    """
    s = wcs_suffix
    if not has_needed_keywords(canvas_hdu.header, s):
        raise ValueError(f"canvas_hdu must include an appropriate WCS: {s}")

    f = quantity_from_table("flux", table, default_unit=u.Unit("ph s-1"))
    if s == "D":
        x = quantity_from_table("x_mm", table, default_unit=u.mm).to(u.mm)
        y = quantity_from_table("y_mm", table, default_unit=u.mm).to(u.mm)
    else:
        arcsec = u.arcsec
        canvas_unit = u.Unit(canvas_hdu.header[f"CUNIT1{s}"])
        x = quantity_from_table("x", table, default_unit=arcsec.to(canvas_unit))
        y = quantity_from_table("y", table, default_unit=arcsec.to(canvas_unit))
        x *= arcsec.to(canvas_unit)
        y *= arcsec.to(canvas_unit)

    xpix, ypix = val2pix(canvas_hdu.header, x.value, y.value, s)

    naxis1 = canvas_hdu.header["NAXIS1"]
    naxis2 = canvas_hdu.header["NAXIS2"]
    # Occasionally 0 is returned as ~ -1e-11
    eps = -1e-7
    mask = (xpix >= eps) * (xpix < naxis1) * (ypix >= eps) * (ypix < naxis2)

    if sub_pixel:
        canvas_hdu = _add_subpixel_sources_to_canvas(
            canvas_hdu, xpix, ypix, f, mask)
    else:
        canvas_hdu = _add_intpixel_sources_to_canvas(
            canvas_hdu, xpix, ypix, f, mask)

    return canvas_hdu


def _add_intpixel_sources_to_canvas(canvas_hdu, xpix, ypix, flux, mask):
    canvas_hdu.header["comment"] = f"Adding {len(flux)} int-pixel files"
    for xpx, ypx, flx, msk in zip(xpix.astype(int), ypix.astype(int),
                                  flux, mask):
        # To prevent adding array values in this manner.
        assert not isinstance(xpx, Iterable), "xpx should be integer"
        canvas_hdu.data[ypx, xpx] += flx.value * msk

    return canvas_hdu


def _add_subpixel_sources_to_canvas(canvas_hdu, xpix, ypix, flux, mask):
    canvas_hdu.header["comment"] = f"Adding {len(flux)} sub-pixel files"
    canvas_shape = canvas_hdu.data.shape
    for xpx, ypx, flx, msk in zip(xpix, ypix, flux, mask):
        if msk:
            xx, yy, fracs = sub_pixel_fractions(xpx, ypx)
            for x, y, frac in zip(xx, yy, fracs):
                if y < canvas_shape[0] and x < canvas_shape[1]:
                    # To prevent adding array values in this manner.
                    assert not isinstance(x, Iterable), "x should be integer"
                    canvas_hdu.data[y, x] += frac * flx.value

    return canvas_hdu


def sub_pixel_fractions(x, y):
    """
    Make a list of pixel coordinates and weights to reflect sub-pixel shifts.

    A point source which isn't centred on a pixel can be modelled by a centred
    PSF convolved with a shifted delta function. A fraction of the delta
    function in moved into each of the adjoining pixels. For example, a star
    at ``(x,y)=(0.2, 0.2)`` would be represented by a following pixel weights::

       ---------------
       | 0.16 | 0.04 |
       ---------------
       | 0.64 | 0.16 |
       ---------------

    where (0,0) is the centre of the bottom-left pixel

    Given (x,y) pixel coordinates, this function returns the fractions of flux
    that should go into the surrounding pixels, as well as the coordinates of
    those neighbouring pixels.

    Parameters
    ----------
    x, y : float

    Returns
    -------
    x_pix, y_pix, fracs : list of (int, int, float)
       The x and y pixel coordinates and their corresponding flux fraction

    """
    x0, dx = divmod(x, 1)
    y0, dy = divmod(y, 1)

    xi0 = int(x0)
    xi1 = xi0 + bool(dx)
    yi0 = int(y0)
    yi1 = yi0 + bool(dy)

    f00 = (1. - dx) * (1. - dy)
    f01 = (1. - dx) * dy
    f10 = dx * (1 - dy)
    f11 = dx * dy

    x_pix = [xi0, xi1, xi0, xi1]
    y_pix = [yi0, yi0, yi1, yi1]
    fracs = [f00, f10, f01, f11]

    return x_pix, y_pix, fracs


def overlay_image(small_im, big_im, coords, mask=None, sub_pixel=False):
    """
    Overlay small_im on top of big_im at the position specified by coords.

    ``small_im`` will be centred at ``coords``

    Adapted from:
    ``https://stackoverflow.com/questions/14063070/overlay-a-smaller-image-on-a-larger-image-python-opencv``

    """
    # TODO - Add in a catch for sub-pixel shifts
    if sub_pixel:
        raise NotImplementedError

    # FIXME: this would not be necessary if we used WCS instead of manual 2pix
    coords = np.ceil(np.asarray(coords).round(10)).astype(int)
    y, x = coords.astype(int)[::-1] - np.array(small_im.shape[-2:]) // 2

    # Image ranges
    x1, x2 = max(0, x), min(big_im.shape[-1], x + small_im.shape[-1])
    y1, y2 = max(0, y), min(big_im.shape[-2], y + small_im.shape[-2])

    # Overlay ranges
    x1o, x2o = max(0, -x), min(small_im.shape[-1], big_im.shape[-1] - x)
    y1o, y2o = max(0, -y), min(small_im.shape[-2], big_im.shape[-2] - y)

    # Exit if nothing to do
    if y1 >= y2 or x1 >= x2 or y1o >= y2o or x1o >= x2o:
        return big_im

    if small_im.ndim == 2 and big_im.ndim == 2:
        small_im_3 = small_im[None, :, :]
        big_im_3 = big_im[None, :, :]
    elif small_im.ndim == 3 and big_im.ndim == 3:
        small_im_3 = small_im
        big_im_3 = big_im
    else:
        raise ValueError(f"Dimensions mismatch between big_im and small_im: "
                         f"{big_im.ndim} : {small_im.ndim}")

    if mask is None:
        big_im_3[:, y1:y2, x1:x2] = big_im_3[:, y1:y2, x1:x2] + \
                                    small_im_3[:, y1o:y2o, x1o:x2o]
    else:
        mask = mask[None, y1o:y2o, x1o:x2o] * np.ones(small_im_3.shape[-3])
        mask = mask.astype(bool)
        big_im_3[:, y1:y2, x1:x2][mask] = big_im_3[:, y1:y2, x1:x2][mask] + \
                                          small_im_3[:, y1o:y2o, x1o:x2o][mask]

    return big_im


def rescale_imagehdu(imagehdu: fits.ImageHDU, pixel_scale: float | u.Quantity,
                     wcs_suffix: str = "", conserve_flux: bool = True,
                     spline_order: int = 1) -> fits.ImageHDU:
    """
    Scale the .data array by the ratio of pixel_scale [deg] and CDELTn.

    `pixel_scale` should NOT be passed as a Quantity!

    Parameters
    ----------
    imagehdu : fits.ImageHDU
    pixel_scale : float or Quantity
        the desired pixel scale of the scaled ImageHDU, as applicable to the WCS
        identified by `wcs_suffix` (by default " "). If float the units are assumed
        to be the same as CUNITa; if a Quantity, the unit need to be convertible.
    wcs_suffix : str
        identifier of the WCS to use for rescaling, By default, this is " ".

    conserve_flux : bool

    spline_order : int
        [1..5] Order of the spline interpolation used by
        ``scipy.ndimage.rotate``

    Returns
    -------
    imagehdu : fits.ImageHDU

    """
    # Identify the wcs to which pixel_scale refers to and determine the zoom factor
    wcs_suffix = wcs_suffix or " "
    primary_wcs = WCS(imagehdu.header, key=wcs_suffix[0])

    # make sure that units are correct and zoom factor is positive
    # The length of the zoom factor will be determined by imagehdu.data,
    # which might differ from the dimension of primary_wcs. Here, pick
    # the spatial dimensions only.
    pixel_scale = pixel_scale << u.Unit(primary_wcs.wcs.cunit[0])
    zoom = np.abs(primary_wcs.wcs.cdelt[:2] / pixel_scale.value)

    if len(imagehdu.data.shape) == 3:
        zoom = np.append(zoom, [1.])  # wavelength dimension unscaled if present

    logger.debug("zoom factor: %s", zoom)

    if primary_wcs.naxis != imagehdu.data.ndim:
        # FIXME: this happens often - shouldn't WCSs be trimmed down before? (OC)
        logger.warning("imagehdu.data.ndim is %d, but primary_wcs.naxis with "
                       "key %s is %d, both should be equal.",
                       imagehdu.data.ndim, wcs_suffix, primary_wcs.naxis)

    if all(zoom == 1.):
        # Nothing to do
        return imagehdu

    sum_orig = np.sum(imagehdu.data)

    # Perform the rescaling. Axes need to be inverted because python.
    new_im = ndi.zoom(imagehdu.data, zoom[::-1], order=spline_order)

    if conserve_flux:
        new_im = np.nan_to_num(new_im, copy=False)
        sum_new = np.sum(new_im)
        if sum_new != 0:
            new_im *= sum_orig / sum_new

    imagehdu.data = new_im

    # Rescale all WCSs in the header
    wcses = find_all_wcs(imagehdu.header)
    for ww in wcses:
        if ww.naxis != imagehdu.data.ndim:
            logger.warning("imagehdu.data.ndim is %d, but wcs.naxis with key "
                           "%s is %d, both should be equal.",
                           imagehdu.data.ndim, ww.wcs.alt, ww.naxis)

        if any(ctype != "LINEAR" for ctype in ww.wcs.ctype):
            logger.warning("Non-linear WCS rescaled using linear procedure.")

        # Assuming linearity, a given world coordinate is determined by
        #   VAL = CRVAL  + (PIX  - CRPIX ) * CDELT   (old system)
        #       = CRVAL' + (PIX' - CRPIX') * CDELT'  (new system)
        # CDELT is simply transformed by the zoom factor:
        #   CDELT' = CDELT / ZOOM
        # The transformation keeps CRVAL' = CRVAL, hence
        #   CRPIX' = PIX' - (PIX - CRPIX) * ZOOM
        # The relation between PIX' and PIX is linear
        #   PIX' = CONST + ZOOM * PIX
        # The fix point is PIX = PIX' = 1/2, which is the lower/left edge of the field,
        # thus  PIX' = (1 - ZOOM)/2 + ZOOM * PIX
        # This leads to
        #   CRPIX' = 1/2 + (CRPIX - 1/2) * ZOOM
        #
        # The transformation only applies to spatial coordinates, which we assume to be
        # the first two in the WCS.
        ww.wcs.cdelt[:2] /= zoom[:2]
        ww.wcs.crpix[:2] = 0.5 + (ww.wcs.crpix[:2] - 0.5) * zoom[:2]
        #ww.wcs.crpix[:2] = (zoom[:2] + 1) / 2 + (ww.wcs.crpix[:2] - 1) * zoom[:2]
        logger.debug("new crpix %s", ww.wcs.crpix)

        imagehdu.header.update(ww.to_header())

    return imagehdu


def reorient_imagehdu(imagehdu: fits.ImageHDU, wcs_suffix: str = "",
                      conserve_flux: bool = True,
                      spline_order: int = 1) -> fits.ImageHDU:
    """
    Apply an affine transformation to the image, as given in its header.

    Parameters
    ----------
    imagehdu : fits.ImageHDU

    wcs_suffix : str

    conserve_flux : bool

    spline_order : int
        [1..5] Order of the spline interpolation used by
        ``scipy.ndimage.rotate``

    Returns
    -------
    imagehdu : fits.ImageHDU

    """
    s = wcs_suffix

    hdr = imagehdu.header
    pc_keys = ["PC1_1", "PC1_2", "PC2_1", "PC2_2"]
    if all(key+s in hdr for key in pc_keys) and imagehdu.data is not None:
        xscen, yscen = pix2val(hdr, hdr["NAXIS1"] / 2., hdr["NAXIS2"] / 2., s)
        hdr["CRVAL1" + s] = xscen
        hdr["CRVAL2" + s] = yscen

        mat = np.array([[hdr["PC1_1" + s], hdr["PC1_2" + s], 0],
                        [hdr["PC2_1" + s], hdr["PC2_2" + s], 0],
                        [0,                0,                1]])
        if imagehdu.data.ndim == 2:
            mat = mat[:2, :2]

        new_im = affine_map(imagehdu.data, matrix=mat, reshape=True,
                            spline_order=spline_order)

        if conserve_flux:
            new_im = np.nan_to_num(new_im, copy=False)
            new_im *= np.sum(imagehdu.data) / np.sum(new_im)

        imagehdu.data = new_im
        hdr[f"CRPIX1{s}"] = hdr["NAXIS1"] / 2.
        hdr[f"CRPIX2{s}"] = hdr["NAXIS2"] / 2.
        for card in [f"PC1_1{s}", f"PC1_2{s}", f"PC2_1{s}", f"PC2_2{s}"]:
            hdr.remove(card)
        imagehdu.header = hdr

    elif any("PC1_1" in key for key in imagehdu.header):
        logger.warning(("PC Keywords were found, but not used due to "
                        "different wcs_suffix given: %s \n %s"),
                       wcs_suffix, dict(imagehdu.header))

    return imagehdu


def affine_map(input, matrix=None, rotation_angle: float = 0.,
               shear_angle: float = 0., scale_factor=None,
               reshape: bool = True, spline_order: int = 3):
    """
    Apply an affine transformation matrix to an image around its centre.

    Similar functionality to ``scipy.ndimage.rotate`` but for the
    ``affine_transformation`` function

    Either a 2x2 affine transformation matrix can be supplied, or the rotation,
    shear, and scaling values which are the basis of an affine transformation.


    Parameters
    ----------
    input : array_like
        The input array
    matrix : ndarray, optional
        A 2x2 affine transformation matrix
    rotation_angle : float, optional
        [deg] If matrix==None, a rotation matrix is built from this angle
    shear_angle : float, optional
        [deg] If matrix==None, a y-axis shear matrix is built from this angle
    scale_factor : list, array
        [mx, my] If matrix==None, a scaling matrix is built from this list
    reshape : bool, optional
        If True, the array is re-sized to contain the whole transformed image
    spline_order : int, optional
        Default is 3. Spline interpolation order

    Returns
    -------
    output : array-like
        The new mapping of the image

    """
    if matrix is None:
        d2r = np.pi / 180.
        c = np.cos(rotation_angle * d2r) if rotation_angle != 0. else 1.
        s = np.sin(rotation_angle * d2r) if rotation_angle != 0. else 0.
        t = np.tan(shear_angle * d2r) if shear_angle != 0. else 0.
        sx = scale_factor[0] if scale_factor is not None else 1.
        sy = scale_factor[1] if scale_factor is not None else 1.
        mat = (np.array([[c, s], [-s, c]]) @
               np.array([[1, t], [0, 1]]) @
               np.array([[sx, 0], [0, sy]]))
    else:
        mat = np.array(matrix)

    h, w = np.array(input.shape)[-2:]
    if reshape:
        # Find the new corner coordinates using the py3.5+ matmul operator
        x, y = mat[:2, :2] @ np.array([[0, w, w, 0], [0, 0, h, h]])
        out = np.ceil(np.array([y.max() - y.min(),
                                x.max() - x.min()])).astype(int)
    else:
        out = np.array([h, w])

    c_out = 0.5 * out
    c_in = 0.5 * np.array([h, w])
    offset = c_in - c_out.dot(mat[:2, :2])

    if mat.shape[0] == 3:
        offset = np.r_[[0], offset]
        out = np.r_[input.shape[0], out]

    output = ndi.affine_transform(input, np.rot90(mat, 2),
                                  output_shape=out, offset=offset,
                                  order=spline_order)

    return output


def add_imagehdu_to_imagehdu(image_hdu: fits.ImageHDU,
                             canvas_hdu: fits.ImageHDU,
                             spline_order: int = 1,
                             wcs_suffix: str = "",
                             conserve_flux: bool = True) -> fits.ImageHDU:
    """
    Re-project one ``fits.ImageHDU`` onto another ``fits.ImageHDU``.

    ..assumption:: of equal grid coordinate lengths

    Parameters
    ----------
    image_hdu : fits.ImageHDU
        The ``ImageHDU`` which will be reprojected onto `canvas_hdu`

    canvas_hdu : fits.ImageHDU
        The ``ImageHDU`` onto which the image_hdu should be projected.
        This must include a valid WCS

    spline_order : int, optional
        Default is 1. The order of the spline interpolator used by the
        ``scipy.ndimage`` functions

    wcs_suffix : str or WCS
        To determine which WCS to use. "" for sky HDUs and "D" for
        ImagePlane HDUs. Can also be ``astropy.wcs.WCS`` object.

    conserve_flux : bool
        Default is True. Used when zooming and rotating to keep flux constant.

    Returns
    -------
    canvas_hdu : fits.ImageHDU

    """
    # TODO: Add a catch for projecting a large image onto a small canvas
    # FIXME: This can mutate the input HDUs, investigate this and deepcopy
    #        if necessary.
    # TODO: rename wcs_suffix to wcs, ideally always pass WCS in the future...
    if isinstance(wcs_suffix, WCS):
        canvas_wcs = wcs_suffix
    else:
        wcs_suffix = wcs_suffix or " "
        canvas_wcs = WCS(canvas_hdu.header, key=wcs_suffix, naxis=2)

    if isinstance(image_hdu.data, u.Quantity):
        image_hdu.data = image_hdu.data.value

    assert canvas_wcs.wcs.cdelt[0] == canvas_wcs.wcs.cdelt[1], \
        "canvas must have equal pixel scale on both axes"

    canvas_pixel_scale = float(canvas_wcs.wcs.cdelt[0])
    conv_fac = u.Unit(image_hdu.header[f"CUNIT1{wcs_suffix}"].lower()).to(canvas_wcs.wcs.cunit[0])

    new_hdu = rescale_imagehdu(image_hdu, pixel_scale=canvas_pixel_scale / conv_fac,
                               wcs_suffix=canvas_wcs.wcs.alt,
                               spline_order=spline_order,
                               conserve_flux=conserve_flux)
    # TODO: Perhaps add separately formatted WCS logger?
    # logger.debug("fromrescale %s", WCS(new_hdu.header, key=canvas_wcs.wcs.alt))
    new_hdu = reorient_imagehdu(new_hdu,
                                wcs_suffix=canvas_wcs.wcs.alt,
                                spline_order=spline_order,
                                conserve_flux=conserve_flux)

    img_center = np.array([[new_hdu.header["NAXIS1"],
                            new_hdu.header["NAXIS2"]]])
    img_center = (img_center - 1) / 2

    new_wcs = WCS(new_hdu.header, key=canvas_wcs.wcs.alt, naxis=2)
    sky_center = new_wcs.wcs_pix2world(img_center, 0)
    if new_wcs.wcs.cunit[0] == "deg":
        sky_center = _fix_360(sky_center)
    # logger.debug("canvas %s", canvas_wcs)
    # logger.debug("new %s", new_wcs)
    logger.debug("sky %s", sky_center)
    sky_center *= conv_fac
    pix_center = canvas_wcs.wcs_world2pix(sky_center, 0)
    logger.debug("pix %s", pix_center)

    canvas_hdu.data = overlay_image(new_hdu.data, canvas_hdu.data,
                                    coords=pix_center.squeeze())

    return canvas_hdu


def pix2val(header, x, y, wcs_suffix=""):
    """
    Return the real coordinates [deg, mm] for coordinates from a Header WCS.

    Parameters
    ----------
    header : fits.Header
    x, y : float, list, array
    wcs_suffix : str

    Returns
    -------
    a, b : float, array
        [deg, mm] Real coordinates as given by the Header WCS

    """
    s = wcs_suffix

    pckeys = [key + s for key in ["PC1_1", "PC1_2", "PC2_1", "PC2_2"]]
    if all(key in header for key in pckeys):
        pc11, pc12, pc21, pc22 = (header[key] for key in pckeys)
    else:
        pc11, pc12, pc21, pc22 = 1, 0, 0, 1

    if (pc11 * pc22 - pc12 * pc21) != 1.0:
        logger.error("PC matrix det != 1.0")

    da = header["CDELT1"+s]
    db = header["CDELT2"+s]
    x0 = header["CRPIX1"+s] - 1
    y0 = header["CRPIX2"+s] - 1
    a0 = header["CRVAL1"+s]
    b0 = header["CRVAL2"+s]

    a = a0 + da * ((x - x0) * pc11 + (y - y0) * pc12)
    b = b0 + db * ((x - x0) * pc21 + (y - y0) * pc22)

    return a, b


def val2pix(header, a, b, wcs_suffix=""):
    """
    Return the pixel coordinates for real coordinates [deg, mm] from a WCS.

    Parameters
    ----------
    header : fits.Header
    a, b : float, list, array
        [deg, mm]

    Returns
    -------
    x, y : float, array
        [pixel] Pixel coordinates as given by the Header WCS

    """
    s = wcs_suffix

    if isinstance(header, fits.ImageHDU):
        header = header.header

    pckeys = [key + s for key in ["PC1_1", "PC1_2", "PC2_1", "PC2_2"]]
    if all(key in header for key in pckeys):
        pc11, pc12, pc21, pc22 = (header[key] for key in pckeys)
    else:
        pc11, pc12, pc21, pc22 = 1, 0, 0, 1

    if (pc11 * pc22 - pc12 * pc21) != 1.0:
        logger.error("PC matrix det != 1.0")

    da = float(header["CDELT1"+s])
    db = float(header["CDELT2"+s])
    x0 = float(header["CRPIX1"+s]) - 1
    y0 = float(header["CRPIX2"+s]) - 1
    a0 = float(header["CRVAL1"+s])
    b0 = float(header["CRVAL2"+s])

    x = x0 + 1. / da * ((a - a0) * pc11 - (b - b0) * pc21)
    y = y0 + 1. / db * ((a - a0) * pc12 + (b - b0) * pc22)

    return x, y


def calc_footprint(header, wcs_suffix="", new_unit: str = None):
    """
    Return the sky/detector positions [deg/mm] of the corners of a header WCS.

    TODO: The rest of this docstring is outdated, please update!

    The positions returned correspond to the corners of the header's
    image array, in this order::

        (ra, dec) = (0,0), (w, 0), (w, h), (0, h)
        (x, y) = (0,0), (w, 0), (w, h), (0, h)

    where ``w``, ``h`` are equal to NAXIS1 and NAXIS2 from the header.

    Parameters
    ----------
    header : fits.Header

    wcs_suffix : str
        Letter suffix for the WCS keywords, e.g. CDELT1D for image-plane coords

    Returns
    -------
    x, y : arrays of floats
        [deg or mm] x are the coordinates for pixels [0, w, w, 0]
        [deg or mm] y are the coordinates for pixels [0, 0, h, h]

    """
    wcs_suffix = wcs_suffix or " "

    if isinstance(header, fits.ImageHDU):
        logger.warning("Passing a HDU to calc_footprint will be deprecated "
                       "in v1.0. Please pass the header only.")
        header = header.header

    # TODO: maybe celestial instead??
    coords = WCS(header, key=wcs_suffix, naxis=2)
    if header["NAXIS"] == 3:
        xy1 = coords.calc_footprint(center=False, axes=(header["NAXIS1"],
                                                        header["NAXIS2"]))
    else:
        try:
            xy1 = coords.calc_footprint(center=False)
        except AstropyWarning:
            # This is only relevant if Warnings are treated as Errors,
            # otherwise astropy will silently return None by itself.
            # Yes this whole thing should be improved, but not now...
            xy1 = None

    if xy1 is None:
        x_ext = max(header["NAXIS1"] - 1, 0)
        y_ext = max(header["NAXIS2"] - 1, 0)
        xy0 = np.array([[0, 0], [0, y_ext], [x_ext, y_ext], [x_ext, 0]])
        xy1 = coords.wcs_pix2world(xy0, 0)

    if (cunit := coords.wcs.cunit[0]) == "deg":
        xy1 = _fix_360(xy1)

    if new_unit is not None:
        convf = cunit.to(new_unit)
        xy1 *= convf

    return xy1


def calc_table_footprint(table: Table, x_name: str, y_name: str,
                         tbl_unit: str, new_unit: str,
                         padding=None) -> np.ndarray:
    """
    Equivalent to ``calc_footprint()``, but for tables instead of images.

    Parameters
    ----------
    table : astropy.table.Table
        Table containing data.
    x_name : str
        Name of the column in `table` to use as x-coordinates.
    y_name : str
        Name of the column in `table` to use as y-coordinates.
    tbl_unit : str
        Default unit to use for x and y if no units are found in `table`.
    new_unit : str
        Unit to convert x and y to, can be identical to `tbl_unit`.
    padding : astropy.units.Quantity, optional
        Constant value to subtract from minima and add to maxima. If used, must
        be Quantity with same physical type as x and y. If None (default), no
        padding is added.

    Returns
    -------
    extent : (4, 2) array
        Array containing corner points (clockwise from bottom left). Format and
        order are equivalent to the output of
        ``astropy.wcs.WCS.calc_footprint()``.

    """
    if padding is not None:
        padding = padding.to(new_unit).value
    else:
        padding = 0.

    x_convf = unit_from_table(x_name, table, tbl_unit).to(new_unit)
    y_convf = unit_from_table(y_name, table, tbl_unit).to(new_unit)

    x_col = table[x_name] * x_convf
    y_col = table[y_name] * y_convf

    x_min = x_col.min() - padding
    x_max = x_col.max() + padding
    y_min = y_col.min() - padding
    y_max = y_col.max() + padding

    extent = np.array([[x_min, y_min],
                       [x_min, y_max],
                       [x_max, y_max],
                       [x_max, y_min]])

    return extent


def split_header(hdr, chunk_size, wcs_suffix=""):
    """
    Split a header into many smaller parts of the chunk_size.

    Parameters
    ----------
    hdr
    chunk_size
    wcs_suffix

    Returns
    -------
    hdr_list
    """
    # TODO: test that this works
    s = wcs_suffix
    naxis1, naxis2 = hdr["NAXIS1"+s], hdr["NAXIS2"+s]
    x0_pix, y0_pix = hdr["CRPIX1"+s], hdr["CRPIX2"+s]       # pix
    x0_sky, y0_sky = hdr["CRVAL1"+s], hdr["CRVAL2"+s]       # deg
    x_delt, y_delt = hdr["CDELT1"+s], hdr["CDELT2"+s]       # deg / pix

    hdr_list = []
    for x1_pix in range(0, naxis1, chunk_size):
        for y1_pix in range(0, naxis2, chunk_size):
            x1_sky = x0_sky + (x1_pix - x0_pix) * x_delt
            y1_sky = y0_sky + (y1_pix - y0_pix) * y_delt
            x2_sky = x1_sky + x_delt * min(chunk_size, naxis1 - x1_pix)
            y2_sky = y1_sky + y_delt * min(chunk_size, naxis2 - y1_pix)

            hdr_sky = header_from_list_of_xy([x1_sky, x2_sky],
                                             [y1_sky, y2_sky],
                                             pixel_scale=x_delt, wcs_suffix=s)
            hdr_list.append(hdr_sky)

    return hdr_list


def _fix_360(arr):
    """Fix the "full circle overflow" that occurs with deg."""
    if isinstance(arr, u.Quantity):
        arr[arr > 270*u.deg] -= 360*u.deg
        arr[arr <= -90*u.deg] += 360*u.deg
    else:
        arr = np.asarray(arr)
        arr[arr > 270] -= 360
        arr[arr <= -90] += 360
    return arr


def _get_unit_from_headers(*headers, wcs_suffix: str = "") -> str:
    unit = headers[0][f"CUNIT1{wcs_suffix}"].lower()
    assert all(header[f"CUNIT{i}{wcs_suffix}"].lower() == unit
               for header, i in product(headers, range(1, 3))), \
        [(i, header[f"CUNIT{i}{wcs_suffix}"])
         for header, i in product(headers, range(1, 3))]
    return unit


def det_wcs_from_sky_wcs(sky_wcs: WCS,
                         pixel_scale: float,
                         plate_scale: float,
                         naxis=None) -> tuple[WCS, np.ndarray]:
    """
    Create detector WCS from celestial WCS using pixel and plate scales.

    Parameters
    ----------
    sky_wcs : astropy.wcs.WCS
        Celestial WCS.
    pixel_scale : float
        Quantity or float (assumed to be arcsec / pixel).
    plate_scale : float
        Quantity or float (assumed to be arcsec / mm).
    naxis : (int, int), optional
        Shape of the image, usually ``NAXIS1`` and ``NAXIS2``. If the input WCS
        holds this information, the default None will use that. Otherwise not
        providing `naxis` will raise and error.

    Returns
    -------
    det_wcs : astropy.wcs.WCS
        Detector WCS.
    det_naxis : (int, int)
        Shape of the image (``NAXIS1``, ``NAXIS2``).

    """
    # TODO: Using astropy units for now to avoid deg vs. arcsec confusion.
    #       Once Scopesim is consistent there, remove astropy units.
    pixel_scale <<= u.arcsec / u.pixel
    plate_scale <<= u.arcsec / u.mm
    logger.debug("Pixel scale: %s", pixel_scale)
    logger.debug("Plate scale: %s", plate_scale)

    pixel_size = pixel_scale / plate_scale
    # TODO: add check if cunit is consistent along all axes
    cunit = sky_wcs.wcs.cunit[0]
    corners = sky_wcs.calc_footprint(center=False, axes=naxis) * cunit
    logger.debug("WCS sky corners:\n%s", corners)
    corners /= plate_scale
    corners = corners.to(u.mm)
    logger.debug("WCS det corners:\n%s", corners)

    return create_wcs_from_points(corners, pixel_size, "D")


def sky_wcs_from_det_wcs(det_wcs: WCS,
                         pixel_scale: float,
                         plate_scale: float,
                         naxis=None) -> tuple[WCS, np.ndarray]:
    """
    Create celestial WCS from detector WCS using pixel and plate scales.

    Parameters
    ----------
    det_wcs : astropy.wcs.WCS
        Detector WCS.
    pixel_scale : float
        Quantity or float (assumed to be arcsec / pixel).
    plate_scale : float
        Quantity or float (assumed to be arcsec / mm).
    naxis : (int, int), optional
        Shape of the image, usually ``NAXIS1`` and ``NAXIS2``. If the input WCS
        holds this information, the default None will use that. Otherwise not
        providing `naxis` will raise and error.

    Returns
    -------
    sky_wcs : astropy.wcs.WCS
        Celestial WCS.
    sky_naxis : (int, int)
        Shape of the image (``NAXIS1``, ``NAXIS2``).

    """
    # TODO: Using astropy units for now to avoid deg vs. arcsec confusion.
    #       Once Scopesim is consistent there, remove astropy units.
    pixel_scale <<= u.arcsec / u.pixel
    plate_scale <<= u.arcsec / u.mm
    logger.debug("Pixel scale: %s", pixel_scale)
    logger.debug("Plate scale: %s", plate_scale)

    # TODO: add check if cunit is consistent along all axes
    cunit = det_wcs.wcs.cunit[0]
    corners = det_wcs.calc_footprint(center=False, axes=naxis) * cunit
    logger.debug("WCS det corners: %s", corners)
    corners *= plate_scale
    corners = corners.to(u.arcsec)
    logger.debug("WCS sky corners: %s", corners)

    return create_wcs_from_points(corners, pixel_scale)
