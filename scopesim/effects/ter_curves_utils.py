# -*- coding: utf-8 -*-
"""Auxiliary functions for ter_curves.py."""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.io.votable import parse_single_table
from astropy.io import ascii as ioascii
from synphot import SpectralElement, Empirical1D

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

    .. versionchanged:: 0.11.2

       Added ``fill_value=0.`` to synphot return mode to fix extrapolation.

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
        return SpectralElement(Empirical1D, points=wave, lookup_table=trans, fill_value=0.)
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
                               lookup_table=tbl["transmission"],
                               fill_value=0.)
    elif filter_name in FILTER_DEFAULTS:
        filt = download_svo_filter(FILTER_DEFAULTS[filter_name])
    else:
        try:
            filt = download_svo_filter(filter_name)
        except ConnectionError:
            filt = None

    return filt


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
