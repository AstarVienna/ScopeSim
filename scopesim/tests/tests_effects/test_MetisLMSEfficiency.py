"""Tests for MetisLMSEfficiency effect"""
import os
import pytest
from pytest import approx

from scopesim.effects.ter_curves import TERCurve
from scopesim.effects.metis_lms_trace_list import MetisLMSEfficiency
from scopesim import rc

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/METIS_LMS/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]

# pylint: disable=no-self-use, missing-class-docstring
# pylint: disable=missing-function-docstring

@pytest.fixture(name="efficiency", scope="class")
def fixture_efficiency():
    """Instantiate a MetisLMSEfficiency object"""
    return MetisLMSEfficiency(wavelen=4.2, filename="TRACE_LMS.fits")


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
        assert  eff_trans == approx(expected, rel=1e-4)
