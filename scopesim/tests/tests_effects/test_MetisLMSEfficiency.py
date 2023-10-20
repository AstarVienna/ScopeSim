"""Tests for MetisLMSEfficiency effect"""

import pytest
from pytest import approx
from unittest.mock import patch

from scopesim.effects.ter_curves import TERCurve
from scopesim.effects.metis_lms_trace_list import MetisLMSEfficiency


# pylint: disable=no-self-use, missing-class-docstring
# pylint: disable=missing-function-docstring

@pytest.fixture(scope="class")
def patch_mock_path_metis(mock_dir):
    metis_dir = mock_dir / "METIS_LMS"
    with patch("scopesim.rc.__search_path__", [metis_dir]):
        yield


# @pytest.fixture(name="efficiency", scope="class")
# def fixture_efficiency(mock_dir):
#     """Instantiate a MetisLMSEfficiency object"""
#     # metis_dir = mock_dir / "METIS_LMS"
#     # with patch("scopesim.rc.__search_path__", [metis_dir]):
#     return MetisLMSEfficiency(wavelen=4.2, filename="TRACE_LMS.fits")


@pytest.mark.usefixtures("patch_mock_path_metis")
class TestMetisLMSEfficiency:
    def test_initialises_correctly(self):
        eff = MetisLMSEfficiency(wavelen=4.2, filename="TRACE_LMS.fits")
        assert isinstance(eff, MetisLMSEfficiency)
        assert isinstance(eff, TERCurve)

    @pytest.mark.parametrize("wavelen, order", [(3.65, 30), (4.6, 24)])
    def test_has_correct_order(self, wavelen, order):
        eff = MetisLMSEfficiency(wavelen=wavelen, filename="TRACE_LMS.fits")
        assert eff.meta['order'] == order

    @pytest.mark.parametrize("lam0, lam, expected",
                             [(3.65, 3.60, 0.127914),
                              (4.80, 4.83, 0.578098),
                              (4.10, 4.30, 0.023711)])
    def test_gives_correct_throughput(self, lam0, lam, expected):
        eff = MetisLMSEfficiency(wavelen=lam0, filename="TRACE_LMS.fits")
        eff_trans = eff.surface.transmission(lam * 1e4).value
        assert eff_trans == approx(expected, rel=1e-4)
