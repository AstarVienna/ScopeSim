"""Tests for METIS WCU classes"""

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

import pytest
from unittest.mock import patch

import numpy as np
from astropy import units as u

from scopesim import UserCommands
from scopesim.effects.metis_wcu import WCUSource, FPMask
from scopesim.utils import seq

def _patched_cmds(mode="wcu_lms", wavelen=3.9, bin_width=0.0003):
    """Minimal UserCommands object to stand in for missing OpticalTrain"""
    return UserCommands(properties={"!OBS.modes": mode,
                                    "!OBS.wavelen": wavelen,
                                    "!SIM.spectral.spectral_bin_width": bin_width})

@pytest.fixture(scope="class")
def patch_mock_path_basic_instrument(mock_dir):
    basic_instrument_dir = mock_dir / "basic_instrument"
    with patch("scopesim.rc.__search_path__", [basic_instrument_dir]):
        yield

def _patched_cmds_lss(mode="wcu_lss", filtername="J", bin_width=0.002):
    """Minimal UserCommands object for wcu_lss mode"""
    return UserCommands(properties={"!OBS.modes": mode,
                                    "!OBS.filter_name": filtername,
                                    "!INST.filter_file_format": "filters/TC_filter_{}.dat",
                                    "!SIM.spectral.spectral_bin_width": bin_width})

@pytest.fixture(name="bbsource", scope="function")
def fixture_bbsource():
    return WCUSource(current_lamp="bb",
                     bb_temp=1000*u.K,
                     is_temp=300*u.K,
                     wcu_temp=300*u.K,
                     bb_aperture=1.,
                     bb_to_is=None,
                     rho_tube=0.95,
                     rho_is=0.95,
                     emiss_mask=1.,
                     diam_is=250,
                     diam_is_in=25.4,
                     diam_is_out=100.,
                     emiss_bb=0.98,
                     current_fpmask="open",
                     fpmask_angle=0,
                     fpmask_shift=(0, 0),
                     fpmask_filename_format="fp_mask_{}.dat",
                     cmds=_patched_cmds())

@pytest.mark.usefixtures("patch_mock_path_basic_instrument")
class TestWCUSource:
    def test_wcu_initialises_correctly(self, bbsource):
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

    def test_get_wavelength_for_lss(self, bbsource):
        bbsource.cmds.update(properties={"!OBS.modes": "wcu_lss",
                                         "!SIM.spectral.spectral_bin_width": 0.002,
                                         "!OBS.filter_name": "J",
                                         "!INST.filter_file_format": "filters/TC_filter_{}.dat"})
        print(bbsource.cmds["!OBS.modes"])
        bbsource.get_wavelength()
        lam = seq(1.15, 1.37, 0.002)
        assert np.all(bbsource.wavelength.value == lam)

    def test_bb_aperture_initialises_correctly(self, bbsource):
        assert bbsource.bb_aperture == 1.

    def test_bb_aperture_set_good_value(self, bbsource):
        bbsource.set_bb_aperture(0.3)
        assert bbsource.bb_aperture == 0.3

    def test_bb_aperture_set_clip_negative_value(self, bbsource):
        bbsource.set_bb_aperture(-3)
        assert bbsource.bb_aperture == 0

    def test_bb_aperture_set_clip_large_value(self, bbsource):
        bbsource.set_bb_aperture(13)
        assert bbsource.bb_aperture == 1

    @pytest.mark.parametrize("newvalue", [1., 0.6, 0.])
    def test_bb_aperture_changes_lamp_emission(self, bbsource, newvalue):
        bbsource.set_bb_aperture(1.)
        ref_lamp = bbsource.intens_lamp
        bbsource.set_bb_aperture(newvalue)
        assert np.allclose(bbsource.intens_lamp, newvalue * ref_lamp)

    @pytest.mark.parametrize("newvalue", [1., 0.6, 0.])
    def test_bb_aperture_does_not_change_background_emission(self, bbsource,
                                                             newvalue):
        bbsource.set_bb_aperture(1.)
        ref_bg = bbsource.intens_bg
        bbsource.set_bb_aperture(newvalue)
        assert np.all(bbsource.intens_bg == ref_bg)


@pytest.fixture(name="fpmask", scope="function")
def fixture_fpmask(mock_path):
    return FPMask(maskname=str(mock_path / "fp_mask_pinhole.dat"))

@pytest.fixture(name="openmask", scope="function")
def fixture_openmask():
    return FPMask(maskname="open")

@pytest.fixture(name="pinholemask", scope="function")
def fixture_pinholemask(mock_path):
    return FPMask(maskname="pinhole",
                  fpmask_filename_format=str(mock_path / "fp_mask_{}.dat"))

class TestFPMask:
    def test_fpmask_initialises_correctly(self, fpmask):
        assert isinstance(fpmask, FPMask)

    def test_fpmask_open_has_correct_hdus(self, openmask):
        assert openmask.holehdu.data is None
        assert openmask.opaquehdu is None

    @pytest.mark.usefixtures("no_file_error")
    def test_fpmask_uses_file_format(self, fpmask, pinholemask):
        assert fpmask.data_container.meta['filename'] == \
            pinholemask.data_container.meta['filename']
        assert np.all(fpmask.holehdu.data == pinholemask.holehdu.data)
        assert np.all(fpmask.opaquehdu.data == pinholemask.opaquehdu.data)

    def test_has_table(self, fpmask):
        assert fpmask.data_container.table is not None

    def test_has_holehdu(self, fpmask):
        assert fpmask.holehdu is not None

    def test_has_opaquehdu(self, fpmask):
        assert fpmask.opaquehdu is not None

    def test_pixarea_correct(self, fpmask):
        hdr = fpmask.holehdu.header
        pixarea = hdr['CDELT1'] * hdr['CDELT2'] * u.arcsec**2
        assert fpmask.pixarea == pixarea

    def test_data_correct(self, fpmask):
        assert fpmask.holehdu.data[241, 1943] == 0
        assert fpmask.holehdu.data[1023, 1023] < np.pi * (0.007532**2) / 4
        assert fpmask.holehdu.data[1021:1025, 1021:1025].sum() == np.pi * (0.007532**2) / 4
        assert fpmask.opaquehdu.data[1023, 1023] == 0
        assert fpmask.opaquehdu.data[748, 1308] == fpmask.pixarea.value
