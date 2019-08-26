from os import path as pth
import numpy as np

from astropy import units as u
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.io import ascii as ioascii
from synphot import SpectralElement, Empirical1D, SourceSpectrum, ConstFlux1D, \
    Observation
from synphot.units import PHOTLAM

from ..utils import find_file, quantity_from_table
from ..rc import __pkg_dir__


FILTER_DEFAULTS = {"U": "Generic/Bessell.U",
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


def download_svo_filter(filter_name):
    """
    Query the SVO service for the true transmittance for a given filter

    Copied 1 to 1 from tynt by Brett Morris

    Parameters
    ----------
    filter_name : str
        Name of the filter as available on the spanish VO filter service
        e.g: ``Paranal/HAWKI.Ks``

    Returns
    -------
    filt_curve : ``synphot.SpectralElement``
        Astronomical filter object.

    """
    path = download_file('http://svo2.cab.inta-csic.es/'
                         'theory/fps3/fps.php?ID={}'.format(filter_name),
                         cache=True)

    true_transmittance = Table.read(path, format='votable')
    wave = true_transmittance['Wavelength'].data.data * u.Angstrom
    trans = true_transmittance['Transmission'].data.data
    filt_curve = SpectralElement(Empirical1D, points=wave, lookup_table=trans)

    return filt_curve


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
        except:
            filt = None

    return filt


def get_zero_mag_spectrum(system_name="AB"):
    if system_name.lower() in ["vega"]:
        spec = vega_spectrum()
    elif system_name.lower() in ["ab"]:
        spec = ab_spectrum()
    elif system_name.lower() in ["st", "hst"]:
        spec = st_spectrum()

    return spec


def vega_spectrum(mag=0):
    vega = SourceSpectrum.from_file(pth.join(__pkg_dir__, "vega.fits"))
    return vega * 10**(-0.4 * mag)


def ab_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.ABmag)


def st_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.STmag)


def zero_mag_flux(filter_name, photometric_system, return_filter=False):
    """
    Returns the zero magnitude photon flux for a filter

    Acceptable filter names are those given in
    ``scopesim.effects.ter_curves_utils.FILTER_DEFAULTS`` or a string with an
    appropriate name for a filter in the Spanish-VO filter-service. Such strings
    must use the naming convention: observatory/instrument.filter. E.g:
    ``paranal/HAWKI.Ks``, or ``Gemini/GMOS-N.CaT``.

    Parameters
    ----------
    filter_name : str
        Name of the filter - see above

    photometric_system : str
        ["vega", "AB", "ST"] Name of the photometric system

    return_filter : bool, optional
        If True, also returns the filter curve object

    Returns
    -------
    flux : float
        [``PHOTLAM``]
    filt : ``synphot.SpectralElement``
        If ``return_filter`` is True

    """

    filt = get_filter(filter_name)
    spec = get_zero_mag_spectrum(photometric_system)

    obs = Observation(spec, filt)
    flux = obs.effstim(flux_unit=PHOTLAM)

    if return_filter:
        return flux, filt
    else:
        return flux

