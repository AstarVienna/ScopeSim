import pytest
import os
from synphot import SpectralElement, SourceSpectrum

from scopesim.effects import SkycalcTERCurve
from scopesim import rc
from scopesim.utils import from_currsys

if rc.__config__["!SIM.tests.run_skycalc_ter_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring SkyCalc integration tests")

FILES_PATH = os.path.join(os.path.dirname(__file__), "../MOCKS/files/")
if FILES_PATH not in rc.__search_path__:
    rc.__search_path__ += [FILES_PATH]


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

    def test_initialises_with_include_false_has_no_table(self):
        sky_ter = SkycalcTERCurve(include=False)
        assert sky_ter.skycalc_table is None

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

    def test_initialise_with_local_skycalc_file(self):
        sky_ter = SkycalcTERCurve(use_local_skycalc_file="skycalc_override.fits")
        assert sky_ter.skycalc_table is not None
