from os import path as pth

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.utils.decorators import deprecated_renamed_argument

from synphot import SourceSpectrum, ConstFlux1D, Empirical1D
from synphot.units import PHOTLAM

from scopesim.rc import __pkg_dir__
from .source import Source
from .. import rc

__all__ = ["empty_sky", "star", "star_field"]


def empty_sky(flux=0):
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
    sky = Source(lam=[0.3, 3.0]*u.um, spectra=[flux, flux]*PHOTLAM,
                 x=[0], y=[0], ref=[0], weight=[1])
    return sky

@deprecated_renamed_argument('mag', 'flux', '0.1.5')
def star(x=0, y=0, flux=0):
    """
    Source object for a single star in either vega, AB magnitudes, or Jansky

    The star is associated with the reference spectrum for each photometric
    system, therefore a reference wavelength or filter does not need to be given

    Parameters
    ----------
    x, y : float
        [arcsec] position from centre of field of view
    flux : float
        [vega mag, AB mag, Jy] Stellar brightness

    Returns
    -------
    src : Source
        A source object with a single entry table field and a reference spectrum

    """

    mag_unit = u.mag
    spec_template = vega_spectrum
    if isinstance(flux, u.Quantity):
        if flux.unit.physical_type == "spectral flux density":  # ABmag and Jy
            mag_unit = u.ABmag
            spec_template = ab_spectrum
            flux = flux.to(u.ABmag)
        flux = flux.value

    spec = spec_template()

    w = 10**(-0.4 * flux)
    ref = 0

    tbl = Table(data=[[x], [y], [w], [ref], [flux]],
                names=["x", "y", "weight", "ref", "mag"],
                units=[u.arcsec, u.arcsec, None, None, mag_unit])
    tbl.meta["photometric_system"] = "vega" if mag_unit == u.mag else "ab"
    src = Source(spectra=spec, table=tbl)

    return src


def star_field(n, mmin, mmax, width, height=None, use_grid=False):
    """
    Creates a super basic field of stars with random positions and brightnesses

    Parameters
    ----------
    n : int
        number of stars

    mmin, mmax : float, astropy.Quantity
        [mag, ABmag, Jy] min and max magnitudes/fluxes of the population stars.
        If floats, then assumed Quantity is vega magnitudes

    width : float
        [arcsec] width of region to put stars in

    height : float, optional
        [arcsec] if None, then height=width

    use_grid : bool, optional
        Place stars randomly or on a grid

    Returns
    -------
    stars : scopesim.Source object
        A Source object with a field of stars that can be fed into the method:
        OpticalTrain.observe()

    See Also
    --------
    OpticalTrain.observe
    OpticalTrain.readout

    """
    if height is None:
        height = width

    mag_unit = u.mag
    spec_template = vega_spectrum
    if isinstance(mmin, u.Quantity):
        if mmin.unit.physical_type == "spectral flux density":  # ABmag and Jy
            mag_unit = u.ABmag
            spec_template = ab_spectrum
            mmin, mmax = mmin.to(u.ABmag), mmax.to(u.ABmag)
        mmin, mmax = mmin.value, mmax.value

    spec = spec_template()

    if rc.__config__["!SIM.random.seed"] is not None:
        np.random.seed(rc.__config__["!SIM.random.seed"])

    if use_grid:
        nw = np.ceil(n**0.5)
        nh = np.ceil(n / nw)
        x, y = np.mgrid[0:1:1/nw, 0:1:1/nh] - 0.5
        positions = x.flatten(), y.flatten()
    else:
        positions = np.random.random(size=(2, n)) - 0.5

    x = width * positions[0]
    y = height * positions[1]
    mags = np.linspace(mmin, mmax, n)
    w = 10**(-0.4 * mags)
    ref = np.zeros(n, dtype=int)

    tbl = Table(data=[x, y, w, ref, mags],
                names=["x", "y", "weight", "ref", "mag"],
                units=[u.arcsec, u.arcsec, None, None, mag_unit])
    tbl.meta["photometric_system"] = "vega" if mag_unit == u.mag else "ab"
    stars = Source(spectra=spec, table=tbl)

    return stars


def vega_spectrum(mag=0):
    vega = SourceSpectrum.from_file(pth.join(__pkg_dir__, "vega.fits"))
    vega = vega * 10 ** (-0.4 * mag)
    return vega


def st_spectrum(mag=0):
    # ..todo: the waves vector is a bit random, in particular its length, but sets the resolution of
    #         the final spectrum in scopesim. Can this be make more general?
    waves = np.geomspace(100, 300000, 50000)
    sp = ConstFlux1D(amplitude=mag*u.STmag)

    return SourceSpectrum(Empirical1D, points=waves, lookup_table=sp(waves))


def ab_spectrum(mag=0):
    # ..todo: the waves vector is a bit random, in particular its length, but sets the resolution of
    #         the final spectrum in scopesim. Can this be make more general?
    waves = np.geomspace(100, 300000, 50000)
    sp = ConstFlux1D(amplitude=mag * u.ABmag)

    return SourceSpectrum(Empirical1D, points=waves, lookup_table=sp(waves))
