import os
from synphot import SpectralElement, SourceSpectrum

from scopesim.effects import SkycalcTERCurve
from scopesim import rc

TRAVIS = True if "TRAVIS" in os.environ else False


def teardown_module():
    if os.path.exists("skycalc_temp.fits"):
        os.remove("skycalc_temp.fits")


class TestInit:
    def test_initialises_with_nothing(self):
        if TRAVIS:
            assert isinstance(SkycalcTERCurve(), SkycalcTERCurve)
            assert False

    def test_initialises_with_some_kwargs(self):
        if TRAVIS:
            sky_ter = SkycalcTERCurve(pwv=1.0, observatory="paranal")
            assert str(sky_ter.surface.meta["wavelength_unit"]) == "um"

    def test_initialises_with_bang_strings(self):
        if TRAVIS:
            rc.__currsys__["!OBS.pwv"] = 20.0
            sky_ter = SkycalcTERCurve(pwv="!OBS.pwv")
            assert sky_ter.skycalc_conn.values["pwv"] == 20.0
            assert isinstance(sky_ter.surface.transmission, SpectralElement)
            assert isinstance(sky_ter.surface.emission, SourceSpectrum)
