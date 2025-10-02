"""Tests for the AtmoLibraryTERCurve class"""

from unittest.mock import patch
import pytest
from scopesim.effects import AtmoLibraryTERCurve
import logging
LOGGER = logging.getLogger(__name__)

MOCKFILE = "test_AtmoLibraryTERCurve.fits"

# pylint: disable=missing-class-docstring,
# pylint: disable=missing-function-docstring

class TestLocalFile:
    def test_initialises_from_local_file(self, mock_path):
        with patch("scopesim.rc.__search_path__", [mock_path]):
            atmo = AtmoLibraryTERCurve(filename=MOCKFILE,
                                       pwv=1)
            assert isinstance(atmo, AtmoLibraryTERCurve)

    @pytest.mark.parametrize("pwv, extname",
                             [(1.0, "PWV_01"),
                              (23, "PWV_23"),
                              (26.5, "PWV_24"),
                              (99,  "PWV_50")])
    def test_picks_nearest_pwv(self, pwv, extname, mock_path):
        with patch("scopesim.rc.__search_path__", [mock_path]):
            atmo = AtmoLibraryTERCurve(filename=MOCKFILE,
                                       pwv=pwv)
            assert atmo.meta['extname'] == extname
            assert atmo.meta['pwv'] == pwv

    def test_updates_with_new_pwv(self, mock_path):
        with patch("scopesim.rc.__search_path__", [mock_path]):
            oldpwv, newpwv = 1, 22.5
            atmo = AtmoLibraryTERCurve(filename=MOCKFILE,
                                       pwv=oldpwv)
            assert atmo.meta['extname'] == "PWV_01"
            assert atmo.meta['pwv'] == oldpwv
            atmo.update(pwv=newpwv)
            assert atmo.meta['extname'] == "PWV_23"
            assert atmo.meta['pwv'] == newpwv

    def test_update_with_unknown_parameter_issues_warning(self, mock_path, caplog):
        caplog.set_level(logging.WARNING)
        with patch("scopesim.rc.__search_path__", [mock_path]):
            atmo = AtmoLibraryTERCurve(filename=MOCKFILE,
                                       pwv=1.)
            atmo.update(temp=23)
        assert "Can only update with parameter pwv" in caplog.text
