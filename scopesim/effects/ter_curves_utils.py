# -*- coding: utf-8 -*-
"""Auxiliary functions for ter_curves.py."""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.io.votable import parse_single_table
from astropy.io import ascii as ioascii
from astropy.io import fits
from synphot import SpectralElement, SourceSpectrum, Empirical1D, Observation
from synphot.units import PHOTLAM

from scopesim.source.source_templates import vega_spectrum, st_spectrum, \
    ab_spectrum
from ..utils import find_file, quantity_from_table, get_logger
from ..server.download_utils import create_client, handle_download


logger = get_logger(__name__)


FILTER_DEFAULTS = {
    "U": "Generic/Bessell.U",
    "B": "Generic/Bessell.B",
    "V": "Generic/Bessell.V",
    "R": "Generic/Bessell.R",
    "I": "Generic/Bessell.I",
    "J": "2MASS/2MASS.J",
    "H": "2MASS/2MASS.H",
    "Ks": "2MASS/2MASS.Ks",
    "K": "Generic/Johnson_UBVRIJHKL.K",
    "L": "Gemini/NIRI.Lprime-G0207w",
    "M": "Gemini/NIRI.Mprime-G0208w",
    "N": "Generic/Johnson_UBVRIJHKL.N",
    "u": "SLOAN/SDSS.u",
    "g": "SLOAN/SDSS.g",
    "r": "SLOAN/SDSS.r",
    "i": "SLOAN/SDSS.i",
    "z": "SLOAN/SDSS.z",
    "u'": "SLOAN/SDSS.uprime_filter",
    "g'": "SLOAN/SDSS.gprime_filter",
    "r'": "SLOAN/SDSS.rprime_filter",
    "i'": "SLOAN/SDSS.iprime_filter",
    "z'": "SLOAN/SDSS.zprime_filter",
    "HAlpha": "Gemini/GMOS-N.Ha",
    "PaBeta": "Gemini/NIRI.PaBeta-G0221",
    "BrGamma": "Gemini/NIRI.BrG-G0218",
}

PATH_HERE = Path(__file__).parent
PATH_SVO_DATA = PATH_HERE.parent / "data" / "svo"


def get_filter_effective_wavelength(filter_name):
    # TODO: This is technically stored in the SVO XML file as WavelengthEff...
    # (actually WavelengthMean, by definition of formula ...)
    if not isinstance(filter_name, str):
        return filter_name

    assert FILTER_DEFAULTS.get(
        filter_name), f"{filter_name} not found in FILTER_DEFAULTS"
    wave, trans = download_svo_filter(
        FILTER_DEFAULTS[filter_name], return_style="quantity")
    eff_wave = (wave * trans).sum() / trans.sum() << u.um

    return eff_wave


def _handle_svo_download(filename: str, params: dict) -> Path:
    # The SVO is only accessible over http, not over https.
    # noinspection HttpUrlsUsage
    base_url = "http://svo2.cab.inta-csic.es/theory/fps3/"
    path = find_file(filename, path=[PATH_SVO_DATA], silent=True)

    # TODO: Turn this into try-except once error_on_missing_file can be True
    #       by default. Actually, check if there are any other places where
    #       error_on_missing_file applies other than this module...
    if not path:
        logger.debug("File not found in %s, downloading...", PATH_SVO_DATA)
        # TODO: Implement proper caching for non-standard filter files.
        path = PATH_SVO_DATA / filename
        client = create_client(base_url)
        handle_download(client, "fps.php", path, params=params)

    return path

def download_svo_filter(filter_name, return_style="synphot"):
    """
    Query the SVO service for the true transmittance for a given filter.

    Adapted from tynt by Brett Morris.

    Parameters
    ----------
    filter_name : str
        Name of the filter as available on the spanish VO filter service
        e.g: "Paranal/HAWKI.Ks"

    return_style : str, optional
        Defines the format the data is returned
        - "synphot": ``synphot.SpectralElement``
        - "table": ``astropy.table.Table``
        - "quantity": ``astropy.unit.Quantity`` [wave, trans]
        - "array": ``np.ndarray`` [wave, trans], where `wave` is in Angstrom
        - "vo_table": ``astropy.io.votable.tree.Table`` - original output from
        SVO service

    Returns
    -------
    filt_curve : See return_style
        Astronomical filter object.

    """
    path = _handle_svo_download(f"{filter_name}.xml", {"ID": filter_name})

    try:
        # tbl = Table.read(path, format="votable")
        votbl = parse_single_table(path)
    except ValueError:
        logger.error("Unable to load %s from %s.", filter_name, path)
        raise

    if return_style == "vo_table":
        return votbl

    tbl_meta = _parse_votable_params(votbl)
    wave = u.Quantity(votbl.array["Wavelength"].data,
                      tbl_meta["WavelengthUnit"], copy=False)
    trans = votbl.array["Transmission"].data

    if return_style == "synphot":
        return SpectralElement(Empirical1D, points=wave, lookup_table=trans)
    if return_style == "table":
        filt = Table(data=[wave, trans], names=["wavelength", "transmission"])
        filt.meta["wavelength_unit"] = str(wave.unit)
        filt.meta["votable_meta"] = tbl_meta  # Don't pollute actual meta...
        return filt
    if return_style == "quantity":
        return wave, trans
    if return_style == "array":
        return wave.value, trans
    raise ValueError(f"return_style {return_style} unknown.")


def download_svo_filter_list(observatory, instrument, short_names=False,
                             include=None, exclude=None):
    """
    Query the SVO service for a list of filter names for an instrument.

    Parameters
    ----------
    observatory : str
        Name of the observatory as available on the spanish VO filter service
        e.g: "Paranal/HAWKI.Ks" --> Paranal

    instrument : str
        Name of the instrument. Be careful of hyphens etc. e.g. "HAWK-I".

    short_names : bool
        Default False. If True, the full SVO names (obs/inst.filt) are split to
        only return the (filt) part of the name.

    include, exclude: str
        Each a string sequence for excluding or including specific filters
        E.g. GTC/OSIRIS has curves for "sdss_g" and "sdss_g_filter".
        We can force the inclusion of only the filter curves by setting
        ``include="_filter"``.

    Returns
    -------
    names : list
        A list of filter names

    """
    path = _handle_svo_download(
        f"{observatory}/{instrument}.xml",
        {"Facility": observatory, "Instrument": instrument})

    tbl = Table.read(path, format="votable")
    names = list(tbl["filterID"])
    if short_names:
        names = [name.split(".")[-1] for name in names]
    if include is not None:
        names = [name for name in names if include in name]
    if exclude is not None:
        names = [name for name in names if exclude not in name]

    return names


def get_filter(filter_name):
    # first check locally
    # check generics
    # check spanish_vo
    path = find_file(filter_name, silent=True)

    if path is not None:
        tbl = ioascii.read(path)
        wave = quantity_from_table("wavelength", tbl, u.um).to(u.um)
        filt = SpectralElement(Empirical1D, points=wave,
                               lookup_table=tbl["transmission"])
    elif filter_name in FILTER_DEFAULTS:
        filt = download_svo_filter(FILTER_DEFAULTS[filter_name])
    else:
        try:
            filt = download_svo_filter(filter_name)
        except ConnectionError:
            filt = None

    return filt


def get_zero_mag_spectrum(system_name: str = "AB") -> SourceSpectrum:
    """
    Return a synphot spectrum of 0 mag in given photometric system.

    Parameters
    ----------
    system_name : {"AB", "Vega", "ST", "HST"}
        Name of the photometric system.

    Raises
    ------
    ValueError
        Raised if `system_name` is invalid.

    Returns
    -------
    spec : synphot.SourceSpectrum
        Output spectrum.

    """
    if system_name.lower() == "vega":
        return vega_spectrum()
    if system_name.lower() == "ab":
        return ab_spectrum()
    if system_name.lower() in ("st", "hst"):
        return st_spectrum()
    raise ValueError(f"system_name {system_name} is unknown")


def zero_mag_flux(filter_name, photometric_system, return_filter=False):
    """
    Return the zero magnitude photon flux for a filter.

    Acceptable filter names are those given in
    ``scopesim.effects.ter_curves_utils.FILTER_DEFAULTS`` or a string with an
    appropriate name for a filter in the Spanish-VO filter-service. Such
    strings must use the naming convention: observatory/instrument.filter. E.g:
    "paranal/HAWKI.Ks", or "Gemini/GMOS-N.CaT".

    Parameters
    ----------
    filter_name : str
        Name of the filter - see above.

    photometric_system : {"vega", "AB", "ST"}
        Name of the photometric system

    return_filter : bool, optional
        If True, also returns the filter curve object

    Returns
    -------
    flux : float
        [PHOTLAM]
    filt : ``synphot.SpectralElement``
        If `return_filter` is True

    """
    filt = get_filter(filter_name)
    spec = get_zero_mag_spectrum(photometric_system)

    obs = Observation(spec, filt)
    flux = obs.effstim(flux_unit=PHOTLAM)

    if return_filter:
        return flux, filt
    return flux


def scale_spectrum(spectrum, filter_name, amplitude):
    """
    Scale a SourceSpectrum to a value in a filter.

    Parameters
    ----------
    spectrum : synphot.SourceSpectrum

    filter_name : str
        Name of a filter from
        - a local instrument package (available in ``rc.__search_path__``)
        - a generic filter name (see ``ter_curves_utils.FILTER_DEFAULTS``)
        - a spanish-vo filter service reference (e.g. "Paranal/HAWKI.Ks")

    amplitude : astropy.Quantity, float
        The value that the spectrum should have in the given filter. Acceptable
        astropy quantities are:
        - u.mag : Vega magnitudes
        - u.ABmag : AB magnitudes
        - u.STmag : HST magnitudes
        - u.Jy : Jansky per filter bandpass
        Additionally the "FLAM" and "FNU" units from ``synphot.units`` can
        be used when passing the quantity for `amplitude`.

    Returns
    -------
    spectrum : synphot.SourceSpectrum
        Input spectrum scaled to the given amplitude in the given filter.

    Examples
    --------
    ::

        >>> from scopesim.source.source_templates import vega_spectrum
        >>> from scopesim.effects.ter_curves_utils as ter_utils
        >>>
        >>> spec = vega_spectrum()
        >>> vega_185 = ter_utils.scale_spectrum(spec, "Ks", -1.85 * u.mag)
        >>> ab_0 = ter_utils.scale_spectrum(spec, "Ks", 0 * u.ABmag)
        >>> jy_3630 = ter_utils.scale_spectrum(spec, "Ks", 3630 * u.Jy)

    """
    if isinstance(amplitude, u.Quantity):
        if amplitude.unit.physical_type == "spectral flux density":
            if amplitude.unit != u.ABmag:
                amplitude = amplitude.to(u.ABmag)
            ref_spec = ab_spectrum(amplitude.value)

        elif amplitude.unit.physical_type == "spectral flux density wav":
            if amplitude.unit != u.STmag:
                amplitude = amplitude.to(u.STmag)
            ref_spec = st_spectrum(amplitude.value)

        elif amplitude.unit == u.mag:
            ref_spec = vega_spectrum(amplitude.value)

        else:
            raise ValueError(f"Units of amplitude must be one of "
                             f"[u.mag, u.ABmag, u.STmag, u.Jy]: {amplitude}")
    else:
        ref_spec = vega_spectrum(amplitude)

    filt = get_filter(filter_name)
    ref_flux = Observation(ref_spec, filt).effstim(flux_unit=PHOTLAM)

    real_flux = Observation(spectrum, filt).effstim(flux_unit=PHOTLAM)
    scale_factor = ref_flux / real_flux
    spectrum *= scale_factor.value

    return spectrum


def apply_throughput_to_cube(
    cube: fits.ImageHDU,
    thru: SpectralElement | SourceSpectrum,
    wave_cube: u.Quantity,
) -> fits.ImageHDU:
    """
    Apply throughput curve to a spectroscopic cube.

    Parameters
    ----------
    cube : fits.ImageHDU
        Three-dimensional image, dimension 0 (in python convention) is the
        spectral dimension.
    thru : SpectralElement | SourceSpectrum
        Throughput curve, spectrum or filter.
    wave_cube : u.Quantity["length"]
        Wavelength axis of the cube.

    Returns
    -------
    cube : fits.ImageHDU
        Header unchanged, data multiplied with wavelength-dependent throughput.

    """
    # Note: wave_cube used to be converted to AA, but synphot understands um
    #       or whatever just as well...
    cube.data *= thru(wave_cube).value[:, None, None]
    return cube


def combine_two_spectra(spec_a, spec_b, action, wave_min, wave_max):
    """
    Combine transmission and/or emission spectrum with a common waverange.

    Spec_A is the source spectrum
    Spec_B is either the transmission or emission that should be applied

    Parameters
    ----------
    spec_a : synphot.SourceSpectrum
    spec_b : synphot.SpectralElement, synphot.SourceSpectrum
    action: {"multiply", "add"}
    wave_min, wave_max : quantity
        [Angstrom]

    Returns
    -------
    new_source : synphot.SourceSpectrum

    """
    if spec_a.waveset is None:
        wave_val = spec_b.waveset.value
    else:
        wave_val = spec_a.waveset.value
    mask = (wave_val > wave_min.value) * (wave_val < wave_max.value)

    wave = ([wave_min.value] + list(wave_val[mask]) + [wave_max.value]) * u.AA
    if "mult" in action.lower():
        spec_c = spec_a(wave) * spec_b(wave)
        # Diagnostic plots - not for general use
        # from matplotlib import pyplot as plt
        # plt.plot(wave, spec_a(wave), label="spec_a")
        # plt.plot(wave, spec_b(wave), label="spec_b")
        # plt.plot(wave, spec_c, label="spec_c")
        # plt.xlim(2.9e4, 4.2e4)
        # plt.legend()
        # plt.show()
    elif "add" in action.lower():
        spec_c = spec_a(wave) + spec_b(wave)
    else:
        raise ValueError(f"action {action} unknown")

    new_source = SourceSpectrum(Empirical1D, points=wave, lookup_table=spec_c)
    new_source.meta.update(spec_b.meta)
    new_source.meta.update(spec_a.meta)

    return new_source


def add_edge_zeros(tbl, wave_colname):
    if isinstance(tbl, Table):
        vals = np.zeros(len(tbl.colnames))
        col_i = np.where(col == wave_colname for col in tbl.colnames)[0][0]
        sgn = np.sign(np.diff(tbl[wave_colname][:2]))
        vals[col_i] = tbl[wave_colname][0] * (1 - 1e-7 * sgn)
        tbl.insert_row(0, vals)
        vals[col_i] = tbl[wave_colname][-1] * (1 + 1e-7 * sgn)
        tbl.insert_row(len(tbl), vals)

    return tbl


def _parse_votable_params(votable):
    """Convert VO Table XML PARAM fields to dict, use Quantity if possible."""
    def _get_metadict(table):
        for param in votable.params:
            if param.datatype == "char" or param.unit is None:
                yield param.name, param.value
            else:
                yield param.name, u.Quantity(param.value, str(param.unit))
    return dict(_get_metadict(votable))
