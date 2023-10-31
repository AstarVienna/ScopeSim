from pathlib import Path

import pytest
from unittest.mock import patch
from synphot import SpectralElement, SourceSpectrum

from scopesim.effects import SkycalcTERCurve
from scopesim import rc

if rc.__config__["!SIM.tests.run_skycalc_ter_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring SkyCalc integration tests")


# TODO: is this at all needed?
@pytest.fixture(scope="module", autouse=True)
def setup_and_teardown():
    Path(rc.__config__["!SIM.file.local_packages_path"]).mkdir(parents=True,
                                                               exist_ok=True)
    yield
    Path("skycalc_temp.fits").unlink(missing_ok=True)


@pytest.mark.filterwarnings("ignore::astropy.units.core.UnitsWarning")
class TestInit:
    @pytest.mark.webtest
    def test_initialises_with_nothing(self):
        assert isinstance(SkycalcTERCurve(), SkycalcTERCurve)

    def test_initialises_with_include_false_has_no_table(self):
        sky_ter = SkycalcTERCurve(include=False)
        assert sky_ter.skycalc_table is None

    @pytest.mark.webtest
    def test_initialises_with_some_kwargs(self):
        sky_ter = SkycalcTERCurve(pwv=1.0, observatory="paranal")
        assert str(sky_ter.surface.meta["wavelength_unit"]) in ("nm", "um")

    @pytest.mark.webtest
    def test_initialises_with_bang_strings(self):
        patched = {"!OBS.pwv": 20.0}
        with patch.dict("scopesim.rc.__currsys__", patched):
            sky_ter = SkycalcTERCurve(pwv="!OBS.pwv")

        assert sky_ter.skycalc_conn.values["pwv"] == 20.0
        assert isinstance(sky_ter.surface.transmission, SpectralElement)
        assert isinstance(sky_ter.surface.emission, SourceSpectrum)

    @pytest.mark.webtest
    def test_initialises_with_non_skycalc_keys(self):
        sky_ter = SkycalcTERCurve(name="bogus")
        assert "name" not in sky_ter.skycalc_conn.values

    def test_initialise_with_local_skycalc_file(self, mock_path):
        sky_ter = SkycalcTERCurve(
            use_local_skycalc_file=str(mock_path / "skycalc_override.fits"))
        assert sky_ter.skycalc_table is not None
