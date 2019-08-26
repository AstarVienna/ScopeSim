import pytest
from pytest import approx

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from astropy import units as u

import synphot as sp

from scopesim.effects import ter_curves_utils as terutils


def test_all_zero_spectra_line_up():
    mag = 0
    vega = terutils.vega_spectrum(mag)
    ab = terutils.ab_spectrum(mag)
    st = terutils.st_spectrum(mag)

    wave = 0.55 * u.um
    assert st(wave).value == approx(vega(wave).value, rel=0.03)
    assert ab(wave).value == approx(vega(wave).value, rel=0.03)


class TestFunctionGetFilter:
    @pytest.mark.parametrize("filt_name", ["V", "Ks", "L", "z'"])
    def test_returns_generic_filter_from_svo(self, filt_name):
        trans = terutils.get_filter(filt_name)
        wave = np.logspace(-1, 1, 1000) * u.um
        assert np.max(trans(wave)) > 0.9

    def test_returns_specific_from_svo(self):
        ks = terutils.get_filter("Paranal/HAWKI.BrGamma")
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
    new_flux = terutils.zero_mag_flux(filter_name, phot_system)
    assert flux == approx(new_flux.value, rel=0.08)



