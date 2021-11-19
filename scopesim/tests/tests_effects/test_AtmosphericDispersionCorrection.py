import pytest
from pytest import approx
import numpy as np

from scopesim.effects import AtmosphericDispersionCorrection as ADC
from scopesim.effects import AtmosphericDispersion as AD
from scopesim.tests.mocks.py_objects.fov_objects import _centre_fov


@pytest.fixture(scope="function")
def atmo_params():
    """
    airmass = 1.14
    altitude = 2400m
    temperature = 7 deg
    pressure = 0.755 bar
    humidity = ?

    Approx atmospheric refraction at 500nm = 24.8 arcsec
    Diff atmo refr relative to 500nm
    - 0.5um : 0 arcsec
    - 1.5um : -0.49 arcsec
    - 2.5um : -0.53 arcsec
    """
    _atmo_params = {"airmass": 1.14,
                    "temperature": 7,
                    "humidity": 0.5,
                    "pressure": 0.755,
                    "latitude": -26,
                    "altitude": 2400,
                    "pupil_angle": 0,
                    "pixel_scale": 1,
                    "wave_min": 0.5,
                    "wave_mid": 0.5,
                    "wave_max": 2.5}
    return _atmo_params


@pytest.mark.usefixtures("atmo_params")
class TestInit:
    def test_initialises_when_all_needed_keywords_given(self, atmo_params):
        assert isinstance(ADC(**atmo_params), ADC)

    def test_throws_error_when_not_all_keywords_are_provided(self):
        with pytest.raises(ValueError):
            ADC(**{"its_over": 9000})


@pytest.mark.usefixtures("atmo_params")
class TestApplyTo:
    def test_does_nothing_when_passed_wrong_type(self, atmo_params):
        adc = ADC(**atmo_params)
        doc_brown = adc.apply_to({"gigawatts": 1.21})
        assert doc_brown["gigawatts"] == 1.21

    def test_zero_shift_at_zenith(self, atmo_params):
        fov = _centre_fov(n=50, waverange=[1.5, 1.7])
        old_crpix_d = fov.header["CRPIX1D"], fov.header["CRPIX2D"]
        atmo_params["airmass"] = 1.
        adc = ADC(**atmo_params)
        fov_new = adc.apply_to(fov)
        new_crpix_d = fov_new.header["CRPIX1D"], fov_new.header["CRPIX2D"]

        assert np.all(old_crpix_d == new_crpix_d)

    @pytest.mark.parametrize("waves, offset",
                             [([0.4, 0.6], 0.),
                              ([1.4, 1.6], 0.490),
                              ([2.4, 2.6], 0.528)])
    def test_correct_test_shift_applied_to_image_plane_wcs(self, atmo_params,
                                                           waves, offset):
        fov = _centre_fov(n=10, waverange=waves)
        old_crpix_d = fov.header["CRPIX1D"], fov.header["CRPIX2D"]

        adc = ADC(**atmo_params)
        fov_new = adc.apply_to(fov)

        new_crpix_d = fov_new.header["CRPIX1D"], fov_new.header["CRPIX2D"]

        abs_diff = np.sum((np.array(new_crpix_d) -
                           np.array(old_crpix_d))**2)**0.5

        # this works because the pixel_scale is 1 arcsec
        assert abs_diff == approx(offset, rel=1e-3)
        assert new_crpix_d[1] == approx(old_crpix_d[1] - offset, rel=1e-3)


@pytest.mark.usefixtures("atmo_params")
class TestCombinedWithAtmoDisp:
    @pytest.mark.parametrize("waves", [(0.7, 0.8), (1.4, 1.6), (2.4, 2.6)])
    @pytest.mark.parametrize("angle", [0, 15, 45, 85, 90])
    @pytest.mark.parametrize("pixel_scale", [0.004, 0.04, 0.4])
    def test_shifts_between_adc_and_ad_are_opposite(self, atmo_params, waves,
                                                    angle, pixel_scale):
        fov_wave_mid = np.average(waves)
        atmo_params["pixel_scale"] = pixel_scale
        atmo_params["pupil_angle"] = angle
        atmo_params["sub_pixel_fraction"] = 0.001

        fov = _centre_fov(n=10, waverange=waves)
        fov.header["CDELT1"] = 1 / 3600 * pixel_scale
        fov.header["CDELT2"] = 1 / 3600 * pixel_scale
        old_crpix_d = np.array([fov.header["CRPIX1D"], fov.header["CRPIX2D"]])

        ad = AD(**atmo_params)
        adc = ADC(**atmo_params)
        ad_shifts = ad.fov_grid()
        ad_x_shift = np.interp(fov_wave_mid, ad_shifts[0], ad_shifts[1])
        ad_y_shift = np.interp(fov_wave_mid, ad_shifts[0], ad_shifts[2])

        adc.apply_to(fov)
        new_crpix_d = np.array([fov.header["CRPIX1D"], fov.header["CRPIX2D"]])
        fov_shifts = new_crpix_d - old_crpix_d
        adc_x_shift = fov_shifts[0] * fov.header["CDELT1"] * 3600
        adc_y_shift = fov_shifts[1] * fov.header["CDELT1"] * 3600

        assert adc_x_shift == approx(ad_x_shift, rel=1e-3)
        assert adc_y_shift == approx(ad_y_shift, rel=1e-3)


