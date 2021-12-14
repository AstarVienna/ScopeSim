import logging
from copy import deepcopy

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
from synphot import SourceSpectrum, Empirical1D

from scopesim import utils, rc
from scopesim.optics import image_plane_utils as imp_utils


def is_field_in_fov(fov_header, field, wcs_suffix=""):
    """
    Returns True if Source.field footprint is inside the FieldOfView footprint

    Parameters
    ----------
    fov_header : fits.Header
        Header from a FieldOfView object
    field : [astropy.Table, astropy.ImageHDU]
        Field object from a Source object
    wcs_suffix : str
        ["S", "D"] Coordinate system: Sky or Detector

    Returns
    -------
    is_inside_fov : bool

    """

    if isinstance(field, fits.ImageHDU) and \
            field.header.get("BG_SRC") is not None:
        is_inside_fov = True
    else:
        if isinstance(field, Table):
            x = list(utils.quantity_from_table("x", field,
                                               u.arcsec).to(u.deg).value)
            y = list(utils.quantity_from_table("y", field,
                                               u.arcsec).to(u.deg).value)
            s = wcs_suffix
            cdelt = utils.quantify(fov_header["CDELT1" + s], u.deg).value
            field_header = imp_utils.header_from_list_of_xy(x, y, cdelt, s)
        elif isinstance(field, (fits.ImageHDU, fits.PrimaryHDU)):
            field_header = field.header
        else:
            logging.warning("Input was neither Table nor ImageHDU: {}"
                          "".format(field))
            return False

        ext_xsky, ext_ysky = imp_utils.calc_footprint(field_header, wcs_suffix)
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
    spline_order = utils.from_currsys("!SIM.computing.spline_order")
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
                temp_hdu = imp_utils.add_imagehdu_to_imagehdu(
                    field, temp_hdu, spline_order, wcs_suffix)
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
        field_new = extract_area_from_table(field, fov_volume)
    elif isinstance(field, fits.ImageHDU):
        field_new = extract_area_from_imagehdu(field, fov_volume)
    else:
        raise ValueError("field must be either Table or ImageHDU: {}"
                         "".format(type(field)))

    return field_new


def extract_area_from_table(table, fov_volume):
    """
    Extracts the entries of a Table that fits inside the fov_volume

    Parameters
    ----------
    table : fits.ImageHDU
        The field ImageHDU, either an image of a wavelength [um] cube
    fov_volume : dict
        Contains {"xs": [xmin, xmax], "ys": [ymin, ymax],
                  "waves": [wave_min, wave_max],
                  "xy_unit": "deg" or "mm", "wave_unit": "um"}

    Returns
    -------
    new_imagehdu : fits.ImageHDU

    """
    fov_unit = u.Unit(fov_volume["xy_unit"])
    fov_xs = (fov_volume["xs"] * fov_unit).to(table["x"].unit)
    fov_ys = (fov_volume["ys"] * fov_unit).to(table["y"].unit)

    mask = (table["x"].data >= fov_xs[0].value) * \
           (table["x"].data <  fov_xs[1].value) * \
           (table["y"].data >= fov_ys[0].value) * \
           (table["y"].data <  fov_ys[1].value)
    table_new = table[mask]

    return table_new


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
    if x0p == x1p: x1p += 1
    if y0p == y1p: y1p += 1

    new_hdr = imp_utils.header_from_list_of_xy([x0s, x1s], [y0s, y1s],
                                               pixel_scale=hdr["CDELT1"])

    if hdr["NAXIS"] == 3:

        # Look 0.5*wdel past the fov edges in each direction to catch any
        # slices where the middle wavelength value doesn't fall inside the
        # fov waverange, but up to 50% of the slice is actually inside the
        # fov waverange:
        # E.g. FOV: [1.92, 2.095], HDU bin centres: [1.9, 2.0, 2.1]
        # CDELT3 = 0.1, and HDU bin edges: [1.85, 1.95, 2.05, 2.15]
        # So 1.9 slice needs to be multiplied by 0.3, and 2.1 slice should be
        # multipled by 0.45 to reach the scaled contribution of the edge slices
        # This scaling factor is:
        # f = ((hdu_bin_centre - fov_edge [+/-] 0.5 * cdelt3) % cdelt3) / cdelt3

        hdu_waves = get_cube_waveset(hdr)
        wdel = hdr["CDELT3"]
        wunit = u.Unit(hdr.get("CUNIT3", "AA"))
        fov_waves = utils.quantify(fov_volume["waves"], u.um).to(wunit).value
        mask = ((hdu_waves > fov_waves[0] - 0.5 * wdel) *
                (hdu_waves <= fov_waves[1] + 0.5 * wdel))  # need to go [+/-] half a bin

        # if min(hdu_waves) > min(fov_waves) or max(hdu_waves) < max(fov_waves):
        #     raise ValueError(f"FOV waveset is not a subset of cube waveset: "
        #                      f"{fov_waves} --> {hdu_waves}")

        # OC [2021-12-14] if fov range is not covered by the source return nothing
        if not np.any(mask):
            print("FOV {} um - {} um: not covered by Source".format(fov_waves[0], fov_waves[1]))
            return None

        i0p, i1p = np.where(mask)[0][0], np.where(mask)[0][-1]
        f0 = (abs(hdu_waves[i0p] - fov_waves[0] + 0.5 * wdel) % wdel) / wdel    # blue edge
        f1 = (abs(hdu_waves[i1p] - fov_waves[1] - 0.5 * wdel) % wdel) / wdel    # red edge
        data = imagehdu.data[i0p:i1p+1, y0p:y1p, x0p:x1p]
        data[0, :, :] *= f0
        if i1p > i0p:
            data[-1, :, :] *= f1

        # w0, w1 : the closest cube wavelengths outside the fov edge wavelengths
        # fov_waves : the fov edge wavelengths
        # f0, f1 : the scaling factors for the blue and red edge cube slices
        #
        # w0, w1 = hdu_waves[i0p], hdu_waves[i1p]
        # print(f"\nw0: {w0}, f0: {f0}, {fov_waves}, f1: {f1}, w1: {w1}")

        new_hdr.update({"NAXIS": 3,
                        "NAXIS3": data.shape[0],
                        "CRVAL3": hdu_waves[i0p],
                        "CRPIX3": 0,
                        "CDELT3": hdr["CDELT3"],
                        "CUNIT3": hdr["CUNIT3"],
                        "CTYPE3": hdr["CTYPE3"],
                        "BUNIT":  hdr["BUNIT"]})

    else:
        data = imagehdu.data[y0p:y1p, x0p:x1p]
        new_hdr["SPEC_REF"] = hdr.get("SPEC_REF")

    new_imagehdu = fits.ImageHDU(data=data)
    new_imagehdu.header.update(new_hdr)

    return new_imagehdu


def get_cube_waveset(hdr, return_quantity=False):
    wval, wdel, wpix, wlen, = [hdr[kw] for kw in ["CRVAL3", "CDELT3",
                                                  "CRPIX3", "NAXIS3"]]
    # ASSUMPTION - cube wavelength is in regularly spaced units of um
    wmin = wval - wdel * wpix
    wmax = wmin + wdel * (wlen - 1)
    wunit = u.Unit(hdr.get("CUNIT3", "AA"))

    if "LOG" in hdr.get("CTYPE3", ""):
        hdu_waves = np.logspace(wmin, wmax, wlen)
    else:
        hdu_waves = np.linspace(wmin, wmax, wlen)

    if return_quantity:
        hdu_waves = hdu_waves << wunit
        hdu_waves.to(u.um)

    return hdu_waves


def extract_range_from_spectrum(spectrum, waverange):
    if not isinstance(spectrum, SourceSpectrum):
        raise ValueError(f"spectrum must be of type synphot.SourceSpectrum: "
                         f"{type(spectrum)}")

    wave_min, wave_max = utils.quantify(waverange, u.um).to(u.AA).value
    spec_waveset = spectrum.waveset.to(u.AA).value
    mask = (spec_waveset > wave_min) * (spec_waveset < wave_max)

    if sum(mask) == 0:
        logging.warning(f"Waverange does not overlap with Spectrum waveset: "
                      f"{[wave_min, wave_max]} <> {spec_waveset} "
                      f"for spectrum {spectrum}")
    if wave_min < min(spec_waveset) or wave_max > max(spec_waveset):
        logging.warning(f"Waverange only partially overlaps with Spectrum waveset: "
                      f"{[wave_min, wave_max]} <> {spec_waveset} "
                      f"for spectrum {spectrum}")

    wave = np.r_[wave_min, spec_waveset[mask], wave_max]
    flux = spectrum(wave)

    new_spectrum = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)
    new_spectrum.meta.update(spectrum.meta)

    return new_spectrum


def make_cube_from_table(table, spectra, waveset, fov_header, sub_pixel=False):
    """

    Parameters
    ----------
    table: astropy.Table
    spectra: dict
    waveset: np.ndarray
    fov_header: fits.Header
    sub_pixel: bool, optional

    Returns
    -------
    cube: fits.ImageHDU
        Units of ph/s/m2/bin --> should this be ph/s/m2/um?

    """
    cube = np.zeros((fov_header["NAXIS2"], fov_header["NAXIS1"], len(waveset)))
    dwave = 0.5 * (np.r_[np.diff(waveset), 0] + np.r_[0, np.diff(waveset)])
    # ..todo: dwave is questionable here. What should the FOV cube units be?

    spec_dict = {i: spec(waveset) * dwave for i, spec in spectra.items()}

    cdelt1, cdelt2 = fov_header["CDELT1"], fov_header["CDELT2"]
    crval1, crval2 = fov_header["CRVAL1"], fov_header["CRVAL2"]
    crpix1, crpix2 = fov_header["CRPIX1"], fov_header["CRPIX2"]
    cunit1, cunit2 = fov_header["CUNIT1"], fov_header["CUNIT2"]

    xps = (table["x"].to(cunit1).value - crval1) / cdelt1 + crpix1
    yps = (table["y"].to(cunit2).value - crval2) / cdelt2 + crpix2
    refs, weights = table["ref"], table["weight"]

    for xp, yp, ref, weight in zip(xps, yps, refs, weights):
        cube[int(round(yp)), int(round(xp)), :] += weight * spec_dict[ref].value

    cube = np.rollaxis(cube, 2)

    cdelt3 = np.diff(waveset[:2]).to(u.um)[0]
    hdu = fits.ImageHDU(data=cube)
    hdu.header.update(fov_header)
    hdu.header.update({"CRVAL3": waveset[0].value,
                       "CRPIX3": 0,
                       "CDELT3": cdelt3.value,
                       "CUNIT3": str(cdelt3.unit),
                       "CTYPE3": "WAVE"})

    return hdu
