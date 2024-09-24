# -*- coding: utf-8 -*-
from typing import Optional, Union
from collections.abc import Iterable
from pathlib import Path

import numpy as np
from astropy import wcs, units as u
from astropy.io import fits
from astropy.table import Table
from synphot import SourceSpectrum, Empirical1D, SpectralElement

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


# FIXME: typing: This should work with the more general (and true, since
#        conversion is done anyway) Quantity["length"] (or "area" resp.), but
#        doing so currently causes a NameError. Not sure what's going on.
def photons_in_range(
        spectra: SourceSpectrum,
        wave_min: u.Quantity[u.um] | float,
        wave_max: u.Quantity[u.um] | float,
        area: Optional[Union[u.Quantity[u.m**2], float]] = None,
        bandpass: Optional[SpectralElement] = None,
) -> Union[u.Quantity[u.ph * u.s**-1 * u.m**-2], u.Quantity[u.ph * u.s**-1]]:
    """
    Integrate photons from spectrum in given wavelength range.

    Parameters
    ----------
    spectra : SourceSpectrum
        Input spectrum.
    wave_min : Union[u.Quantity["length"], float]
        Minimum wavelength. If float, assumes um.
    wave_max : Union[u.Quantity["length"], float]
        Maximum wavelength. If float, assumes um.
    area : Optional[Union[u.Quantity["area"], float]], optional
        Area to multiply with. If float, assumes m**2. The default is None.
    bandpass : Optional[SpectralElement], optional
        Filter to take into account, if any. The default is None.

    Returns
    -------
    counts : astropy.units.Quantity
        Either in ph/s/m**2 or just ph/s (if area was given).

    """
    # Note: Assuming um if given as float.
    wave_min = (wave_min << u.um << u.Angstrom).value
    wave_max = (wave_max << u.um << u.Angstrom).value
    # Note: There appear to be some float shenanigans going on here, but
    # rounding produces an error in the spectrum evaluation. Not sure what's
    # going on, maybe it's fine as-is.

    counts = []
    for spec in spectra:
        waveset = spec.waveset.value
        mask = (waveset > wave_min) * (waveset < wave_max)
        wave = np.array([wave_min, *waveset[mask], wave_max])
        flux = spec(wave).value

        # flux [ph s-1 cm-2] == flux [ph s-1 cm-2 AA-1] * wave [AA]
        if isinstance(bandpass, SpectralElement):
            bandpass.model.bounds_error = True
            counts.append(np.trapz(bandpass(wave).value * flux, wave))
        else:
            counts.append(np.trapz(flux, wave))

    # counts = flux [ph s-1 cm-2]
    counts = (counts * u.ph * u.s**-1 * u.cm**-2).to(u.ph * u.s**-1 * u.m**-2)
    if area is not None:
        counts *= (area << u.m**2)

    return counts


def scale_imagehdu(imagehdu, waverange, area=None):
    """
    ..todo: implement this

    For the moment, all imagehdu must be accompanied by a spectrum in PHOTLAM

    Future functionality will include scaling here of:
    ph s-1
    ph s-1 m-2
    ph s-1 m-2
    ph s-1 m-2 um-1
    ph s-1 m-2 um-1 arcsec-2
    J s-1 m-2 Hz-1
    J s-1 m-2 Hz-1 arcsec-2
    ABMAG
    ABMAG arcsec-2
    VEGAMAG
    VEGAMAG arcsec-2
    """
    if "SPEC_REF" not in imagehdu.header:
        raise ValueError("For this version, an ImageHDU must be associated "
                         "with a spectrum. This will change in the future.")

    return imagehdu


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
