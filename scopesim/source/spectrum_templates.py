from os import path as pth

from astropy import units as u
from scopesim.rc import __pkg_dir__
from synphot import SourceSpectrum, ConstFlux1D


def vega_spectrum(mag=0):
    vega = SourceSpectrum.from_file(pth.join(__pkg_dir__, "vega.fits"))
    vega = vega * 10 ** (-0.4 * mag)
    return vega


def st_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.STmag)


def ab_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.ABmag)