"""Tests for METIS WCU classes"""

import pytest
from astropy import units as u

from scopesim.effects.metis_wcu import BlackBodySource
from scopesim.utils import seq

@pytest.fixture(name="bbsource", scope="function")
def fixture_bbsource():
    return BlackBodySource(bb_temp=1000*u.K, wcu_temp=300*u.K)

class TestBlackBodySource:
    def test_initialises_correctly(self, bbsource):
        assert isinstance(bbsource, BlackBodySource)

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
                u.ph / (u.s * u.sr * u.um * u.m**2))

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

    def test_ignore_negative_bb_temperature(self, bbsource):
        old_temp = bbsource.meta['bb_temp']
        bbsource.set_temperature(bb_temp=-1000 * u.K)
        assert bbsource.meta['bb_temp'] == old_temp

    def test_emission_increases_with_temperature(self, bbsource):
        old_temp = bbsource.meta['bb_temp']
        new_temp = 1.2 * old_temp
        old_emission = bbsource.surface.emission
        bbsource.set_temperature(bb_temp=new_temp)
        new_emission = bbsource.surface.emission
        lam_ref = seq(2.2, 15, 0.1) * u.um
        assert all(new_emission(lam_ref) > old_emission(lam_ref))
