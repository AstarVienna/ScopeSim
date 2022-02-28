import pytest
import scopesim.source.source_templates
from pytest import approx

import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u

import scopesim.source.source_templates as src_ts
from scopesim.effects import ter_curves_utils as ter_utils

PLOTS = False


class TestFunctionGetFilter:
    @pytest.mark.parametrize("filt_name", ["V", "Ks", "L", "z'"])
    def test_returns_generic_filter_from_svo(self, filt_name):
        trans = ter_utils.get_filter(filt_name)
        wave = np.logspace(-1, 1, 1000) * u.um
        assert np.max(trans(wave)) > 0.9

    def test_returns_specific_from_svo(self):
        ks = ter_utils.get_filter("Paranal/HAWKI.BrGamma")
        wave = np.linspace(2.1, 2.2, 100) * u.um
        assert np.max(ks(wave)) > 0.75


@pytest.mark.parametrize("filter_name, phot_system, flux",
                         [("V", "vega", 995),
                          ("z", "AB", 602),
                          ("Ks", "vega", 43)])
def test_zero_mag_flux(filter_name, phot_system, flux):
    """
    Photon fluxes are in line with:
    http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    """
    new_flux = ter_utils.zero_mag_flux(filter_name, phot_system)
    assert flux == approx(new_flux.value, rel=0.08)


def test_compare_br_gamma():
    filt = "Ks"
    flux_vega = ter_utils.zero_mag_flux(filt, "vega").value
    flux_ab = ter_utils.zero_mag_flux(filt, "AB").value
    delta_mag = 2.5 * np.log10(flux_ab/flux_vega)
    assert delta_mag == approx(1.85, rel=0.01)


class TestScaleSpectrum:
    def test_scales_vega_spectrum_to_vega_ab_or_jansky(self):
        spec = scopesim.source.source_templates.vega_spectrum()
        vega_185 = ter_utils.scale_spectrum(spec, "Ks", -1.85 * u.mag)
        ab_0 = ter_utils.scale_spectrum(spec, "Ks", 0 * u.ABmag)
        jy_3630 = ter_utils.scale_spectrum(spec, "Ks", 3630 * u.Jy)

        wave = np.linspace(1.8, 2.5, 1000) * u.um
        assert vega_185(wave).value == approx(ab_0(wave).value, rel=1e-2)
        assert vega_185(wave).value == approx(jy_3630(wave).value, rel=1e-2)

        if PLOTS:
            plt.plot(wave, spec(wave), "b")
            plt.plot(wave, vega_185(wave), "r:")
            plt.plot(wave, ab_0(wave), "g--")
            plt.plot(wave, jy_3630(wave), "y-.")
            plt.semilogy()
            plt.show()


class TestGetFilterEffectiveWavelength:
    def test_Ks_is_around_2_2um(self):
        wave_eff = ter_utils.get_filter_effective_wavelength("Ks")
        print(wave_eff)
        assert wave_eff.value == approx(2.19, rel=1e-3)

    def test_rprime_is_around_0_626um(self):
        wave_eff = ter_utils.get_filter_effective_wavelength("r'")
        print(wave_eff)
        assert wave_eff.value == approx(0.626, rel=1e-3)
