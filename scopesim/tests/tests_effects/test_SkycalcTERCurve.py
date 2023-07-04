import pytest
import os
from synphot import SpectralElement, SourceSpectrum
from astropy import units as u

from scopesim.effects import SkycalcTERCurve
from scopesim import rc
from scopesim.utils import save_unit_to

if rc.__config__["!SIM.tests.run_skycalc_ter_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring SkyCalc integration tests")

FILES_PATH = os.path.join(os.path.dirname(__file__), "../mocks/files/")
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


class TestRescalingAndUnitConversion:
    """Test rescaling of SkyCalcTERCurves and unit conversion."""

    def test_convert_to_nm(self):
        """Can we idempotent convert to nm?"""
        tc = SkycalcTERCurve(
            observatory="armazones",
            wmin=0.7,
            wmax=2.5,
            wunit="um",
        )
        tc.convert_to_nm()
        wmin1 = tc.meta["wmin"]
        wmax1 = tc.meta["wmax"]
        wunit1 = tc.meta["wunit"]
        assert wmin1 == 700
        assert wmax1 == 2500
        assert wunit1 == "nm"
        tc.convert_to_nm()
        assert wmin1 == 700
        assert wmax1 == 2500
        assert wunit1 == "nm"

    def test_save_to(self):
        """Test whether exact unit conversion works."""
        for ufrom, uto, fac in [
            (u.um, u.nm, 1000.),
            (u.mm, u.um, 1000.),
            (u.nm, u.um, 0.001),
            (u.um, u.mm, 0.001),
        ]:
            assert save_unit_to(ufrom, uto) == fac

    def test_rescaling(self):
        """Test rescaling of SkyCalcTERCurves and unit conversion.

        There is a unit conversion problem that prevents rescaling of
        SkyCalcTERCurves.

        When wunit="um" and wmax=2.5, then this gets converted to nm by astropy,
        because that is what Skycalc expects. However,
        u.um.to(u.nm) = 999.9999999999999
        and therefore 2.5 um is converted to 2499.99 nm.

        2499.99 nm is less than 2.5 um, so there is no overlap anymore between
        the source spectrum and the to-be-scaled to spectrum, leading to an
        error.
        """
        tc = SkycalcTERCurve(
            observatory="armazones",
            wmin=0.7,
            # It is essential that wmax=2.5 exactly for this test to work properly.
            wmax=2.5,
            wunit="um",
            rescale_emission={
                # The Q filter also goes to 2.5 exactly
                "filter_name": "Q",
                "filename_format": "TC_filter_{}.dat",
                "value": 17,
                "unit": "mag",
            }
        )
        # Accessing tc.emission will lead to an exception:
        # E               synphot.exceptions.PartialOverlap: Source spectrum
        # and bandpass do not fully overlap. You may use force=[extrap|taper]
        # to force this Observation anyway.
        _ = tc.emission
