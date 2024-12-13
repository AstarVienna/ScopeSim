"""Tests for METIS WCU classes"""

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

import pytest
import numpy as np
from astropy import units as u

from scopesim import UserCommands
from scopesim.effects.metis_wcu import WCUSource
from scopesim.utils import seq

def _patched_cmds(mode="wcu_lms", wavelen=3.9, bin_width=0.0003):
    """Minimal UserCommands object to stand in for missing OpticalTrain"""
    return UserCommands(properties={"!OBS.modes": mode,
                                    "!OBS.wavelen": wavelen,
                                    "!SIM.spectral.spectral_bin_width": bin_width})

@pytest.fixture(name="bbsource", scope="function")
def fixture_bbsource():
    return WCUSource(current_lamp="bb",
                     bb_temp=1000*u.K,
                     is_temp=300*u.K,
                     wcu_temp=300*u.K,
                     bb_to_is=None,
                     rho_tube=0.95,
                     rho_is=0.95,
                     rho_mask=0.95,
                     diam_is=250,
                     diam_is_in=25.4,
                     diam_is_out=100.,
                     emiss_bb=0.98,
                     cmds=_patched_cmds())

class TestWCUSource:
    def test_initialises_correctly(self, bbsource):
        assert isinstance(bbsource, WCUSource)

    def test_bbsource_has_temperatures(self, bbsource):
        assert bbsource.meta['bb_temp'] == 1000 * u.K
        assert bbsource.meta['wcu_temp'] == 300 * u.K

    def test_bbsource_has_table(self, bbsource):
        assert bbsource.surface.table

    def test_bbsource_table_has_correct_columns(self, bbsource):
        assert (bbsource.surface.table.colnames ==
                ['wavelength', 'transmission', 'emission'])

    def test_bbsource_has_correct_emission_units(self, bbsource):
        assert (bbsource.surface.table['emission'].unit ==
                u.ph / (u.s * u.arcsec**2 * u.um * u.m**2))

    def test_can_set_bb_temperature(self, bbsource):
        old_temp = bbsource.meta['bb_temp']
        new_temp = 1.5 * old_temp
        bbsource.set_temperature(new_temp)
        assert bbsource.meta['bb_temp'] == new_temp

    def test_can_set_wcu_temperature(self, bbsource):
        old_temp = bbsource.meta['wcu_temp']
        new_temp = old_temp + 12 * u.K
        bbsource.set_temperature(wcu_temp=new_temp)
        assert bbsource.meta['wcu_temp'] == new_temp

    def test_can_set_temperature_in_celsius(self, bbsource):
        bbsource.set_temperature(bb_temp=30*u.deg_C)
        assert bbsource.meta['bb_temp'] == 303.15 * u.K

    def test_ignore_incompatible_units_bb(self, bbsource):
        old_temp = bbsource.meta['bb_temp']
        with pytest.raises(u.UnitConversionError):
            bbsource.set_temperature(bb_temp=2*u.Jy)
        assert bbsource.meta['bb_temp'] == old_temp

    def test_ignore_incompatible_units_wcu(self, bbsource):
        old_temp = bbsource.meta['wcu_temp']
        with pytest.raises(u.UnitConversionError):
            bbsource.set_temperature(wcu_temp=2*u.Tesla)
        assert bbsource.meta['wcu_temp'] == old_temp

    def test_ignore_negative_bb_temperature(self, bbsource):
        old_temp = bbsource.meta['bb_temp']
        with pytest.raises(ValueError):
            bbsource.set_temperature(bb_temp=-1000 * u.K)
        assert bbsource.meta['bb_temp'] == old_temp

    def test_ignore_negative_wcu_temperature(self, bbsource):
        old_temp = bbsource.meta['wcu_temp']
        with pytest.raises(ValueError):
            bbsource.set_temperature(wcu_temp=-1000 * u.K)
        assert bbsource.meta['wcu_temp'] == old_temp

    def test_emission_increases_with_temperature(self, bbsource):
        """
        Test that increasing the temperature leads to increasing emission.

        This is not a rigorous quantitative test, but should stay correct
        when all the wavelength-dependent corrections are made to the
        original black-body emission.
        """
        old_temp = bbsource.meta['bb_temp']
        new_temp = 1.2 * old_temp
        old_emission = bbsource.surface.emission
        bbsource.set_temperature(bb_temp=new_temp)
        new_emission = bbsource.surface.emission
        lam_ref = seq(2.2, 15, 0.1) * u.um
        assert all(new_emission(lam_ref) > old_emission(lam_ref))


    def test_black_body_source_fails_without_parameters(self):
        with pytest.raises(ValueError):
            bbsource = WCUSource()
            assert isinstance(bbsource, WCUSource)


    def test_get_wavelength_for_lms(self, bbsource):
        binw=0.0001
        lamc=3.9
        bbsource.cmds = _patched_cmds(mode="wcu_lms",  wavelen=lamc, bin_width=binw)
        bbsource.get_wavelength()
        lam = 3.6 + (np.arange(6001) * binw)
        assert len(bbsource.wavelength) == len(lam)
        assert np.all(bbsource.wavelength.value == lam)
