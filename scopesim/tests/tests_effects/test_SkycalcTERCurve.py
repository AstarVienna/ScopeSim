import pytest
import os
from synphot import SpectralElement, SourceSpectrum

from scopesim.effects import SkycalcTERCurve
from scopesim import rc
from scopesim.utils import from_currsys

if rc.__config__["!SIM.tests.run_skycalc_ter_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring SkyCalc integration tests")


def setup_module():
    rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
    if not os.path.exists(rc_local_path):
        os.mkdir(rc_local_path)
        rc.__config__["!SIM.file.local_packages_path"] = os.path.abspath(
            rc_local_path)


def teardown_module():
    if os.path.exists("skycalc_temp.fits"):
        os.remove("skycalc_temp.fits")


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(SkycalcTERCurve(), SkycalcTERCurve)

    def test_initialises_with_some_kwargs(self):
        sky_ter = SkycalcTERCurve(pwv=1.0, observatory="paranal")
        assert str(sky_ter.surface.meta["wavelength_unit"]) == "um"

    def test_initialises_with_bang_strings(self):
        rc.__currsys__["!OBS.pwv"] = 20.0
        sky_ter = SkycalcTERCurve(pwv="!OBS.pwv")
        assert sky_ter.skycalc_conn.values["pwv"] == 20.0
        assert isinstance(sky_ter.surface.transmission, SpectralElement)
        assert isinstance(sky_ter.surface.emission, SourceSpectrum)

    def test_initialises_with_non_skycalc_keys(self):
        sky_ter = SkycalcTERCurve(name="bogus")
        assert "name" not in sky_ter.skycalc_conn.values
