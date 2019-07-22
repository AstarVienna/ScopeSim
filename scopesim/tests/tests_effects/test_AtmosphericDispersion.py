import pytest

from scopesim.effects.shifts import AtmosphericDispersion

from ..mocks.py_objects.yaml_objects import _atmo_yaml_dict


@pytest.fixture(scope="function")
def atmo_yaml_dict():
    return _atmo_yaml_dict()


@pytest.mark.usefixtures("atmo_yaml_dict")
class TestInit:
    def test_throws_error_when_initialised_with_nothing(self):
        with pytest.raises(ValueError):
            isinstance(AtmosphericDispersion(), AtmosphericDispersion)

    def test_passes_if_all_keywords_are_present(self, atmo_yaml_dict):
        kwargs_dict = atmo_yaml_dict["properties"]
        atmo_disp = AtmosphericDispersion(**kwargs_dict)

        assert isinstance(atmo_disp, AtmosphericDispersion)

    def test_throws_error_if_one_keyword_missing(self, atmo_yaml_dict):
        kwargs_dict = atmo_yaml_dict["properties"]
        kwargs_dict.pop("airmass")
        with pytest.raises(ValueError):
            isinstance(AtmosphericDispersion(), AtmosphericDispersion)


@pytest.mark.usefixtures("atmo_yaml_dict")
class TestFovGrid:
    def test_returns_similar_values_to_lasilla_website(self, atmo_yaml_dict):
        """
        http: // gtc - phase2.gtc.iac.es / science / astroweb / atmosRefraction.php
        airmass = 1.14
        altitude = 2400m
        temperature = 7 deg
        pressure = 0.755 bar
        humidity = ?

        Approx atmospheric refraction at 500nm = 24.8 arcsec
        Diff atmo refr relative to 500nm
        - 0.5 : 0 arcsec
        - 1.5 : -0.49 arcsec
        - 2.5 : -0.53 arcsec

        """

        kwargs_dict = atmo_yaml_dict["properties"]

        from scopesim import rc
        rc.__config__["!SIM.spectral.lam_min"] = 0.5
        rc.__config__["!SIM.spectral.lam_mid"] = 1.5
        rc.__config__["!SIM.spectral.lam_max"] = 2.5
        kwargs_dict["airmass"] = 1.14
        kwargs_dict["temperature"] = 7
        kwargs_dict["pressure"] = 0.755
        kwargs_dict["altitude"] = 2400
        kwargs_dict["pupil_angle"] = 0

        atmo_disp = AtmosphericDispersion(**kwargs_dict)

        print(atmo_disp.fov_grid())
        # assert atmo_disp.fov_grid()[1] == 0
        # assert atmo_disp.fov_grid()[1] == -0.49
        # assert atmo_disp.fov_grid()[1] == -0.53

        # ..todo:: Fix this too!
