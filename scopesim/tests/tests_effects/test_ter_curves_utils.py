# -*- coding: utf-8 -*-

import pytest
from pytest import approx

import numpy as np

from astropy import units as u

from scopesim.effects import ter_curves_utils as ter_utils

PLOTS = False


# get_filter relies on None from find_file
@pytest.mark.webtest
@pytest.mark.usefixtures("no_file_error")
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


@pytest.mark.webtest
@pytest.mark.usefixtures("no_file_error")
class TestGetFilterEffectiveWavelength:
    def test_Ks_is_around_2_2um(self):
        wave_eff = ter_utils.get_filter_effective_wavelength("Ks")
        print(wave_eff)
        assert wave_eff.value == approx(2.19, rel=1e-3)

    def test_rprime_is_around_0_626um(self):
        wave_eff = ter_utils.get_filter_effective_wavelength("r'")
        print(wave_eff)
        assert wave_eff.value == approx(0.626, rel=1e-3)
