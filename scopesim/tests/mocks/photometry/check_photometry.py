import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
from scopesim.effects import AtmosphericTERCurve, FilterCurve
from scopesim.source.source_templates import empty_sky


def skycalc_plus_eso_filter_gives_similar_photons_flux_to_ohio():
    """
    ohio: http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    skycalc: https://www.eso.org/observing/etc/doc/skycalc/helpskycalc.html#mags

    - Use the skycalc background spectrum
    - Multiply it by a filter curve
    - Sum up the photons
    - Divide by FWHM
    - Rescale to zero mag using skycalc backgorund magnitude

    """

    paranal_bg_mags = {"U": 20.76, "B": 21.14, "V": 20.67, "R": 20.32,
                       "I": 19.48, "Z": 18.67, "Y": 17.58, "J": 16.87,
                       "H": 14.43, "K": 15.23, "L":  6.00, "M":  1.14,
                       "N": -2.29, "Q": -7.27}
    ohio_zero_ph = {"U": 756, "B": 1393, "V": 996, "R": 702, "I": 452,
                    "J": 193, "H": 93, "K": 43.6}       # ph s-1 cm-2 angstrom-1
    ohio_fwhm = {"U": 0.06, "B": 0.09, "V": 0.085, "R": 0.15, "I": 0.15,
                 "J": 0.26, "H": 0.29, "K": 0.41}       # micron

    for filt_name in ohio_zero_ph:

        # load the sky background curve and the filter
        print(filt_name)

        skycalc_file = "eso_skycalc_r_6000_wave_0.3_15nm_airmass_1_pwv_2.5.dat"
        filter_file = f"eso_etc_filter_{filt_name}.dat"
        atmo_ter = AtmosphericTERCurve(filename=skycalc_file)
        filt_ter = FilterCurve(filename=filter_file)

        print(filt_ter.fwhm, filt_ter.centre)

        # get skycalc bg flux
        src = empty_sky()
        src = atmo_ter.apply_to(src)
        src = filt_ter.apply_to(src)

        wave = np.arange(0.2, 2.5, 0.0001) * u.um
        flux = src.spectra[1](wave)

        # report normalised flux ratio in ph s-1 cm-2 angstrom-1 arcsec-2
        ph = np.trapz(flux, wave) / filt_ter.fwhm.to(u.AA)
        ph = ph.to("ph s-1 cm-2 angstrom-1")
        ph_zero_mag = ph / 10 ** (-0.4 * paranal_bg_mags[filt_name])
        ohio_ph = ohio_zero_ph[filt_name]
        print("Sky BG:", ph)
        print(f"Zero mag normalised flux [{ph_zero_mag.unit}]:",
              "Ohio:", ohio_zero_ph[filt_name],
              "ESO:", ph_zero_mag.value)
        print("Ratio:", ph_zero_mag.value / ohio_ph)

        # report absolute flux ratio in ph s-1 cm-2 arcsec-2 (per filter bandpass)
        ph = np.trapz(flux, wave)
        ph = ph.to("ph s-1 cm-2")
        ph_zero_mag = ph / 10 ** (-0.4 * paranal_bg_mags[filt_name])
        ohio_ph = ohio_zero_ph[filt_name] * ohio_fwhm[filt_name] * 1e4
        print("Sky BG:", ph)
        print(f"Zero mag absolute flux [{ph_zero_mag.unit}]:",
              "Ohio:", ohio_ph,
              "ESO:", ph_zero_mag.value)
        print("Ratio:", ph_zero_mag.value / ohio_ph)


def convert_skycalc_fits_to_ascii(filename):
    data = fits.getdata(filename)
    tbl = Table(data)
    tbl.meta["comments"] = ["lam_unit: nm", "flux_unit: ph/s/m2/micron/arcsec2"]
    tbl.write(filename.replace(".fits", ".dat"), format="ascii.fixed_width")

# convert_skycalc_fits_to_ascii("eso_skycalc_r_6000_wave_0.3_15nm_airmass_1_pwv_2.5.fits")

# skycalc_plus_eso_filter_gives_similar_photons_flux_to_ohio()


