from collections.abc import Iterable

import numpy as np
from astropy import wcs, units as u
from astropy.io import fits
from astropy.table import Table
from synphot import SourceSpectrum, Empirical1D, SpectralElement

from ..utils import find_file, quantify, get_logger


logger = get_logger(__name__)


def validate_source_input(**kwargs):
    if "filename" in kwargs and kwargs["filename"] is not None:
        filename = kwargs["filename"]
        if find_file(filename) is None:
            logger.warning("filename was not found: %s", filename)

    if "image" in kwargs and kwargs["image"] is not None:
        image_hdu = kwargs["image"]
        if not isinstance(image_hdu, (fits.PrimaryHDU, fits.ImageHDU)):
            raise ValueError("image must be fits.HDU object with a WCS."
                             f"{type(image_hdu) = }")

        if len(wcs.find_all_wcs(image_hdu.header)) == 0:
            logger.warning("image does not contain valid WCS. %s",
                            wcs.WCS(image_hdu))

    if "table" in kwargs and kwargs["table"] is not None:
        tbl = kwargs["table"]
        if not isinstance(tbl, Table):
            raise ValueError("table must be an astropy.Table object:"
                             f"{type(tbl) = }")

        if not np.all([col in tbl.colnames for col in ["x", "y", "ref"]]):
            raise ValueError("table must contain at least column names: "
                             f"'x, y, ref': {tbl.colnames}")

    return True


def convert_to_list_of_spectra(spectra, lam):
    spectra_list = []
    if isinstance(spectra, SourceSpectrum):
        spectra_list += [spectra]

    elif lam is None and\
            isinstance(spectra, (tuple, list)) and \
            isinstance(spectra[0], SourceSpectrum):
        spectra_list += spectra

    elif lam is not None and len(spectra.shape) == 1 and \
            isinstance(spectra, np.ndarray) and \
            isinstance(lam, np.ndarray):
        spec = SourceSpectrum(Empirical1D, points=lam, lookup_table=spectra)
        spectra_list += [spec]

    elif ((isinstance(spectra, np.ndarray) and
           len(spectra.shape) == 2) or
          (isinstance(spectra, (list, tuple)) and
           isinstance(spectra[0], np.ndarray))) and \
            isinstance(lam, np.ndarray):

        for sp in spectra:
            spec = SourceSpectrum(Empirical1D, points=lam, lookup_table=sp)
            spectra_list += [spec]

    return spectra_list


def photons_in_range(spectra, wave_min, wave_max, area=None, bandpass=None):
    """

    Parameters
    ----------
    spectra
    wave_min
        [um]
    wave_max
        [um]
    area : Quantity
        [m2]
    bandpass : SpectralElement


    Returns
    -------
    counts : u.Quantity array

    """

    if isinstance(wave_min, u.Quantity):
        wave_min = wave_min.to(u.Angstrom).value
    else:
        wave_min *= 1E4

    if isinstance(wave_max, u.Quantity):
        wave_max = wave_max.to(u.Angstrom).value
    else:
        wave_max *= 1E4

    counts = []
    for spec in spectra:
        waveset = spec.waveset.value
        mask = (waveset > wave_min) * (waveset < wave_max)
        x = waveset[mask]
        x = np.append(np.append(wave_min, x), wave_max)
        y = spec(x).value

        # flux [ph s-1 cm-2] == y [ph s-1 cm-2 AA-1] * x [AA]
        if isinstance(bandpass, SpectralElement):
            bp = bandpass(x)
            bandpass.model.bounds_error = True
            counts += [np.trapz(bp * y, x)]
        else:
            counts += [np.trapz(y, x)]

    # counts = flux [ph s-1 cm-2]
    counts = 1E4 * np.array(counts)    # to get from cm-2 to m-2
    counts *= u.ph * u.s**-1 * u.m**-2
    if area is not None:
        counts *= quantify(area, u.m ** 2)

    return counts


def scale_imagehdu(imagehdu, waverange, area=None):
    # ..todo: implement this
    # For the moment, all imagehdu must be accompanied by a spectrum in PHOTLAM
    #
    # Future functionality will include scaling here of:
    # ph s-1
    # ph s-1 m-2
    # ph s-1 m-2
    # ph s-1 m-2 um-1
    # ph s-1 m-2 um-1 arcsec-2
    # J s-1 m-2 Hz-1
    # J s-1 m-2 Hz-1 arcsec-2
    # ABMAG
    # ABMAG arcsec-2
    # VEGAMAG
    # VEGAMAG arcsec-2

    if "SPEC_REF" not in imagehdu.header:
        raise ValueError("For this version, an ImageHDU must be associated "
                         "with a spectrum. This will change in the future.")

    return imagehdu


def make_img_wcs_header(pixel_scale, image_size):
    """
    Create a WCS header for an image

    pixel_scale : float
        arcsecs
    image_size : tuple
        x, y where x, y are integers

    """
    ra, dec = 0, 0
    x, y = image_size

    imgwcs = wcs.WCS(naxis=2)
    imgwcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    imgwcs.wcs.cunit = [u.deg, u.deg]
    imgwcs.wcs.crpix = [(x + 1) / 2, (y + 1) / 2]
    imgwcs.wcs.cdelt = np.array([-pixel_scale / 3600, pixel_scale / 3600])
    imgwcs.wcs.crval = [ra, dec]
    imgwcs.wcs.cunit = [u.deg, u.deg]

    return imgwcs.to_header()
