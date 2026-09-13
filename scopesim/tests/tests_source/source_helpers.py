# -*- coding: utf-8 -*-
"""Aux function for source testing."""

import numpy as np
from scipy import integrate
from astropy import units as u
from synphot import SourceSpectrum, SpectralElement


# This used to be in source_utils, but was ultimately only used in a few tests,
# so moved here as a simple helper function.
def photons_in_range(
    spectra: SourceSpectrum,
    wave_min: u.Quantity[u.um] | float,
    wave_max: u.Quantity[u.um] | float,
    area: u.Quantity[u.m**2] | float | None = None,
    bandpass: SpectralElement | None = None,
) -> u.Quantity[u.ph * u.s**-1 * u.m**-2] | u.Quantity[u.ph * u.s**-1]:
    """
    Integrate photons from spectrum in given wavelength range.

    Parameters
    ----------
    spectra : SourceSpectrum
        Input spectrum.
    wave_min : u.Quantity["length"] or float
        Minimum wavelength. If float, assumes um.
    wave_max : u.Quantity["length"] or float
        Maximum wavelength. If float, assumes um.
    area : u.Quantity["area"] or float, optional
        Area to multiply with. If float, assumes m**2. The default is None.
    bandpass : SpectralElement, optional
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
            counts.append(integrate.trapezoid(bandpass(wave).value * flux, wave))
        else:
            counts.append(integrate.trapezoid(flux, wave))

    # counts = flux [ph s-1 cm-2]
    counts = (counts * u.ph * u.s**-1 * u.cm**-2).to(u.ph * u.s**-1 * u.m**-2)
    if area is not None:
        counts *= (area << u.m**2)

    return counts


# This used to be a method of Source, but was ultimately only used in two tests,
# so moved here as a simple helper function.
def src_photons_in_range(src, wave_min, wave_max, area=None):
    """

    Parameters
    ----------
    wave_min : float, u.Quantity
        [um]
    wave_max : float, u.Quantity
        [um]
    area : float, u.Quantity, optional
        [m2]

    Returns
    -------
    counts : u.Quantity list
        [ph / s / m2] if area is None
        [ph / s] if area is passed

    """
    indices = src.spectra.keys()
    spectra = [src.spectra[ii] for ii in indices]
    counts = photons_in_range(spectra, wave_min, wave_max, area=area,
                              bandpass=src.bandpass)
    return counts
