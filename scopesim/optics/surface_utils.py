
import numpy as np
from astropy import units as u
from synphot import SourceSpectrum, BlackBody1D, Empirical1D

from ..utils import quantify, get_logger


logger = get_logger(__name__)


def make_emission_from_emissivity(temp: u.Quantity[u.K],
                                  emiss_src_spec) -> SourceSpectrum:
    """
    Create an emission SourceSpectrum using blackbody and emissivity curves.

    Parameters
    ----------
    temp : Quantity[Kelvin]
        Blackbody temperature.
    emiss_src_spec : synphot.SpectralElement
        An emissivity response curve in the range [0..1]

    Returns
    -------
    flux : synphot.SourceSpectrum

    """
    if emiss_src_spec is None:
        logger.warning("Either emission or emissivity must be set")
        return None

    # This line is redundant for all the places we actually call this function.
    # But the tests want this to work, so I'll include it. Ultimately, this is
    # an internal utils function, so it should be fine to just to static type
    # checking to ensure temp is always in K.
    with u.set_enabled_equivalencies(u.temperature()):
        temp <<= u.K

    flux = SourceSpectrum(BlackBody1D, temperature=temp.value)
    flux *= emiss_src_spec
    flux.meta["temperature"] = temp
    flux.meta["solid_angle"] = u.sr**-1
    flux.meta["history"] = [
        "Created from Blackbody curve. Units are per steradian",
    ]

    return flux


def make_emission_from_array(flux, wave, meta) -> SourceSpectrum:
    """
    Create an emission SourceSpectrum using an array.

    Takes care of bins and solid angles. The solid_angle is kept in the
    returned SourceSpectrum meta dictionary under self.meta["solid_angle"].

    Parameters
    ----------
    flux : array-like, Quantity
        if flux is not an array, the ``emission_unit`` must be in meta dict
    wave : array-like, Quantity
        if flux is not an array, the ``wavelength_unit`` must be in meta dict
    meta : dict

    Returns
    -------
    flux : synphot.SourceSpectrum

    """
    if not isinstance(flux, u.Quantity):
        try:
            flux = quantify(flux, meta["emission_unit"])
        except KeyError:
            logger.warning("emission_unit must be set in self.meta, "
                           "or emission must be an astropy.Quantity object")
            return None

    if not isinstance(wave, u.Quantity):
        logger.warning("wavelength and emission must be "
                       "astropy.Quantity objects")
        return None

    flux_unit, angle = extract_type_from_unit(flux.unit, "solid angle")
    flux /= angle

    flux = normalise_flux_if_binned(flux, wave)

    orig_unit = flux.unit
    flux = SourceSpectrum(Empirical1D, points=wave,
                          lookup_table=flux)
    flux.meta["solid_angle"] = angle
    flux.meta["history"] = [
        f"Created from emission array with units {orig_unit}",
    ]

    return flux


def normalise_flux_if_binned(flux, wave):
    """
    Convert a binned flux Quantity array back into flux density.

    The flux density normalising unit is taken from the wavelength Quantity
    unit.

    Parameters
    ----------
    flux : array-like Quantity
        flux unit must include ``bin``, e.g. ``ph s-1 m-2 bin-1``
    wave : array-like Quantity

    Returns
    -------
    flux : array-like Quantity

    """
    if (u.bin not in flux.unit.bases and
            "flux density" in str(flux.unit.physical_type)):
        # not binned, return as-is
        return flux

    bins = np.zeros_like(wave)
    # edge bins only have half the flux of other bins
    bins[:-1] = 0.5 * np.diff(wave)
    bins[1:] += 0.5 * np.diff(wave)

    bin_unit = extract_base_from_unit(flux.unit, u.bin)[1]
    # TODO: Why not just do flux /= (bins / u.bin)
    flux /= (bins * bin_unit)

    return flux


# moved these two here from utils, because they weren't used anywhere else
def extract_type_from_unit(unit, unit_type):
    """
    Extract ``astropy`` physical type from a compound unit.

    Parameters
    ----------
    unit : astropy.Unit
    unit_type : str
        The physical type of the unit as given by ``astropy``

    Returns
    -------
    new_unit : Unit
        The input unit minus any base units corresponding to `unit_type`.
    extracted_units : Unit
        Any base units corresponding to `unit_type`.

    """
    extracted_units = u.Unit("")
    for base, power in zip(unit.bases, unit.powers):
        if unit_type == (base ** abs(power)).physical_type:
            extracted_units *= base ** power

    new_unit = unit / extracted_units

    return new_unit, extracted_units


def extract_base_from_unit(unit, base_unit):
    """
    Extract ``astropy`` base unit from a compound unit.

    Parameters
    ----------
    unit : astropy.Unit
    base_unit : Unit, str

    Returns
    -------
    new_unit : Unit
        The input unit minus any base units corresponding to `base_unit`.
    extracted_units : Unit
        Any base units corresponding to `base_unit`.

    """
    extracted_units = u.Unit("")
    for base, power in zip(unit.bases, unit.powers):
        if base == base_unit:
            extracted_units *= base ** power

    new_unit = unit * extracted_units ** -1

    return new_unit, extracted_units
