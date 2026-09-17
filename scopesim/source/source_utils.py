# -*- coding: utf-8 -*-

from collections.abc import Iterable
from pathlib import Path

import numpy as np
from astropy import wcs, units as u
from astropy.io import fits
from astropy.table import Table
from synphot import SourceSpectrum, Empirical1D

from ..utils import find_file, get_logger, convert_table_comments_to_dict


logger = get_logger(__name__)


def validate_source_input(**kwargs) -> None:
    """
    Check validity of kwargs passed to ``Source`` object.

    Currently checks for "filename", "image" and "table", raising the
    exceptions listed below. Additionally logs a warning if no WCS is found in
    an image, or if a given filename cannot be found.

    Parameters
    ----------
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    TypeError
        Raised if an image isn't a FITS HDU or a table isn't an astropy Table.
    ValueError
        Raised if a table does not contain the minimum required columns.

    Returns
    -------
    None

    """
    if (filename := kwargs.get("filename")) is not None:
        if find_file(filename) is None:
            logger.warning("filename was not found: %s", filename)

    if (image_hdu := kwargs.get("image")) is not None:
        if not isinstance(image_hdu, (fits.PrimaryHDU, fits.ImageHDU)):
            raise TypeError(
                f"Image must be fits.HDU object: {type(image_hdu) = }")

        if not wcs.find_all_wcs(image_hdu.header):
            logger.warning(
                "Image does not contain valid WCS. %s", wcs.WCS(image_hdu))

    if (tbl := kwargs.get("table")) is not None:
        if not isinstance(tbl, Table):
            raise TypeError(
                f"Table must be astropy.Table object: {type(tbl) = }")

        if not {"x", "y", "ref"}.issubset(tbl.colnames):
            raise ValueError(
                "Table must contain at least the following column names: 'x', "
                f"""'y', 'ref'; found only: '{"', '".join(tbl.colnames)}'""")
            # TODO py312: The triple quotes will become redundant in 3.12 !


def convert_to_list_of_spectra(spectra, lam) -> list[SourceSpectrum]:
    """Produce SourceSpectrum instances or pass them through."""
    def _synphotify(spec):
        if not isinstance(lam, np.ndarray):
            raise TypeError("If spectra is/are given as array(s), lam must be "
                            "an array as well.")
        return SourceSpectrum(Empirical1D, points=lam, lookup_table=spec)

    def _from_arrays(specarrays):
        for spec in specarrays:
            yield _synphotify(spec)

    def _get_list():
        if isinstance(spectra, SourceSpectrum):
            yield spectra
            return

        if (isinstance(spectra, Iterable) and
                not isinstance(spectra, np.ndarray)):
            _spectra = list(spectra)  # avoid eating iterators in all()
            if all(isinstance(spec, SourceSpectrum) for spec in _spectra):
                yield from _spectra
            elif all(isinstance(spec, np.ndarray) for spec in _spectra):
                yield from _from_arrays(_spectra)
            else:
                raise ValueError(
                    "If given as an iterable, spectra must consist of all "
                    "synphot spectra or all arrays")
            return

        if isinstance(spectra, np.ndarray):
            if spectra.ndim == 1:
                yield _synphotify(spectra)
            elif spectra.ndim == 2:
                yield from _from_arrays(spectra)
            else:
                raise ValueError(
                    "If given as an array, spectra must have either 1 (single "
                    "flux list) or 2 (flux of multiple spectra) dimensions, "
                    f"but {spectra.ndim} were found.")
            return

    return list(_get_list())


def make_img_wcs_header(
    pixel_scale: float,
    image_size: tuple[int, int],
) -> fits.Header:
    """
    Create a WCS header for an image.

    Parameters
    ----------
    pixel_scale : float
        Pixel scale in arcsecs.
    image_size : tuple[int, int]
        Image size (x, y).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    ra, dec = 0, 0
    x, y = image_size

    imgwcs = wcs.WCS(naxis=2)
    imgwcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    imgwcs.wcs.cunit = [u.deg, u.deg]
    imgwcs.wcs.crpix = [(x + 1) / 2, (y + 1) / 2]
    imgwcs.wcs.cdelt = np.array([-pixel_scale, pixel_scale]) / 3600
    imgwcs.wcs.crval = [ra, dec]
    imgwcs.wcs.cunit = [u.deg, u.deg]

    return imgwcs.to_header()


def parse_sed_table(filename: Path | str) -> Table:
    """
    Parse SED table from example cubes.

    Parameters
    ----------
    filename : Path | str
        Input file path.

    Returns
    -------
    astropy.table.Table
        Parsed table.

    """
    tbl = Table.read(filename, format="ascii")
    tbl.meta.update(convert_table_comments_to_dict(tbl))
    tbl.meta.pop("comments")
    new_names = {}
    for col in tbl.columns:
        cmt = tbl.meta[col.replace("col", "column ")].split("(", maxsplit=1)
        tbl[col].unit = cmt[-1].strip(")")
        new_names[col] = cmt[0].split(";", maxsplit=1)[0].strip()
    # Cannot do a single loop because tbl.columns would get mutated...
    for old_name, new_name in new_names.items():
        tbl[old_name].name = new_name
    return tbl
