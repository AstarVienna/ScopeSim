"""Tests for class ADCWheel"""

import pytest
from unittest.mock import patch

from scopesim.effects import TERCurve, ADCWheel


@pytest.fixture(name="adcwheel", scope="function")
def fixture_adcwheel(mock_path):
    """Instantiate an ADCWheel"""
    with patch("scopesim.rc.__search_path__", [mock_path]):
        return ADCWheel(adc_names=["const_90", "const_10"],
                        filename_format="TER_ADC_{}.dat",
                        current_adc="const_90")


# pylint: disable=missing-class-docstring,
# pylint: disable=missing-function-docstring
class TestADCWheel:
    def test_initialises_correctly(self, adcwheel):
        assert isinstance(adcwheel, ADCWheel)

    def test_has_transmission_curves(self, adcwheel):
        assert len(adcwheel.table) == 2

    def test_current_adc_is_ter_curve(self, adcwheel):
        assert isinstance(adcwheel.current_adc, TERCurve)

    def test_current_adc_has_fov_grid_method(self, adcwheel):
        assert hasattr(adcwheel.current_adc, "fov_grid")

    def test_change_to_known_adc(self, adcwheel):
        adcwheel.change_adc('const_10')
        assert adcwheel.current_adc.meta['name'] == 'const_10'

    def test_change_to_unknown_adc(self, adcwheel):
        with pytest.raises(ValueError):
            adcwheel.change_adc('X')

    def test_reports_current_adc_false(self, mock_path):
        with patch("scopesim.rc.__search_path__", [mock_path]):
            adcwheel = ADCWheel(adc_names=["const_90", "const_10"],
                                filename_format="TER_ADC_{}.dat",
                                current_adc=False)
        assert not adcwheel.current_adc

    def test_changes_to_false(self, adcwheel):
        adcwheel.change_adc(False)
        assert not adcwheel.current_adc
