import pytest
from pytest import approx

import numpy as np

from scopesim import rc
from scopesim.effects.shifts import AtmosphericDispersion, \
                                    get_pixel_border_waves_from_atmo_disp

from ..mocks.py_objects.yaml_objects import _atmo_yaml_dict


@pytest.fixture(scope="function")
def atmo_yaml_dict():
    rc.__config__["!SIM.spectral.wave_min"] = 0.5
    rc.__config__["!SIM.spectral.wave_mid"] = 1.5
    rc.__config__["!SIM.spectral.wave_max"] = 2.5
    rc.__config__["!SIM.sub_pixel.fraction"] = 1
    rc.__config__["!INST.pixel_scale"] = 0.004

    atmo_dict = _atmo_yaml_dict()
    atmo_dict["properties"].update({"pupil_angle": 0, "airmass": 1.14,
                                    "altitude": 2400, "temperature": 7,
                                    "pressure": 0.755,
                                    "pixel_scale": 0.004})      # "!INST.pixel_scale"

    return atmo_dict


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
    def test_returns_list_of_3_arrays_with_correct_which(self, atmo_yaml_dict):
        atmo_disp = AtmosphericDispersion(**atmo_yaml_dict["properties"])
        response = atmo_disp.fov_grid()
        assert len(response) == 3
        assert all([isinstance(x, np.ndarray) for x in response])

    def test_returns_non_with_wrong_which_keyword(self, atmo_yaml_dict):
        atmo_disp = AtmosphericDispersion(**atmo_yaml_dict["properties"])
        response = atmo_disp.fov_grid(which="bogus")
        assert response is None

    def test_returns_similar_values_to_lasilla_website(self, atmo_yaml_dict):
        atmo_disp = AtmosphericDispersion(**atmo_yaml_dict["properties"])
        waves, dx, dy = atmo_disp.fov_grid()
        assert dy[0] - dy[-1] == approx(0.53, rel=1e-2)
        assert all(dx == 0)
        assert waves[0] == 0.5 and waves[-1] == 2.5

    def test_returns_same_results_when_turned_90_degrees(self, atmo_yaml_dict):
        atmo_yaml_dict["properties"]["pupil_angle"] = 90
        atmo_disp = AtmosphericDispersion(**atmo_yaml_dict["properties"])
        waves, dx, dy = atmo_disp.fov_grid()
        assert dx[0] - dx[-1]== approx(0.53, rel=1e-2)
        assert all([y == approx(0) for y in dy])

    def test_returns_same_results_when_turned_30_degrees(self, atmo_yaml_dict):
        atmo_yaml_dict["properties"]["pupil_angle"] = 30
        atmo_disp = AtmosphericDispersion(**atmo_yaml_dict["properties"])
        waves, dx, dy = atmo_disp.fov_grid()
        dr = ((dx[0] - dx[-1])**2 + (dy[0] - dy[-1])**2)**0.5
        assert dr == approx(0.53, rel=1e-2)


class TestGetPixelBorderWavesFromAtmoDisp:
    """
    http://gtc-phase2.gtc.iac.es/science/astroweb/atmosRefraction.php
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
    def test_returns_sensible_data(self):
        atmo_params = {"z0"     : 28.7,               # in deg
                       "temp"   : 7,                # in degC
                       "rel_hum": 100,              # in %
                       "pres"   : 755,              # in mbar
                       "lat"    : -26,              # in deg
                       "h"      : 2400,             # in m
                       "wave_min": 0.5,             # in um
                       "wave_mid": 1.5,
                       "wave_max": 2.5,
                       "pixel_scale": 0.004,        # in arcsec
                       "sub_pixel_fraction": 1,
                       "num_steps": 1000}

        waves, shifts = get_pixel_border_waves_from_atmo_disp(**atmo_params)

        assert shifts[0] - shifts[-1] == approx(0.53, rel=2e-2)
        assert waves[0] == 0.5 and waves[-1] == 2.5

        atmo_params["wave_max"] = 1.5
        waves, shifts = get_pixel_border_waves_from_atmo_disp(**atmo_params)
        assert shifts[0] - shifts[-1] == approx(0.49, rel=2e-2)
        assert waves[0] == 0.5 and waves[-1] == 1.5
