"""Tests for class SlitWheel"""

import pytest
from unittest.mock import patch

from scopesim.effects import ApertureMask, SlitWheel


@pytest.fixture(name="swheel", scope="function")
def fixture_swheel(mock_path):
    """Instantiate a SlitWheel"""
    with patch("scopesim.rc.__search_path__", [mock_path]):
        return SlitWheel(slit_names=["A", "B"],
                         filename_format="MASK_slit_{}.dat",
                         current_slit="B")


# pylint: disable=no-self-use, missing-class-docstring,
# pylint: disable=missing-function-docstring
class TestSlitWheel:
    def test_initialises_correctly(self, swheel):
        assert isinstance(swheel, SlitWheel)

    def test_has_aperture_masks(self, swheel):
        assert len(swheel.table) == 2

    def test_current_slit_is_aperture_mask(self, swheel):
        assert isinstance(swheel.current_slit, ApertureMask)

    def test_current_slit_has_fov_grid_method(self, swheel):
        assert hasattr(swheel.current_slit, "fov_grid")

    def test_change_to_known_slit(self, swheel):
        swheel.change_slit('A')
        assert swheel.current_slit.meta['name'] == 'A'

    def test_change_to_unknown_slit(self, swheel):
        with pytest.raises(ValueError):
            swheel.change_slit('X')

    def test_reports_current_slit_false(self, mock_path):
        with patch("scopesim.rc.__search_path__", [mock_path]):
            swheel = SlitWheel(slit_names=["A", "B"],
                               filename_format="MASK_slit_{}.dat",
                               current_slit=False)
        assert not swheel.current_slit

    def test_changes_to_false(self, swheel):
        swheel.change_slit(False)
        assert not swheel.current_slit

    def test_add_slit_to_wheel(self, swheel):
        num_slit_old = len(swheel.slits)
        kwargs = {"array_dict": {"x": [-2, -1, 1, 2],
                                 "y": [-1, -2, 2, 1]},
                  "x_unit": "arcsec",
                  "y_unit": "arcsec"}
        newslit = ApertureMask(name="newslit", **kwargs)

        swheel.add_slit(newslit, name='newslit')
        assert len(swheel.slits) == num_slit_old + 1

        swheel.change_slit('newslit')
        assert swheel.current_slit.display_name == "newslit"
