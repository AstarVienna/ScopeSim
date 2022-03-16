
import os
import pytest

import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u

from scopesim.effects import ter_curves as tc
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim import rc

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]

PLOTS = False

# pylint: disable=no-self-use, missing-class-docstring
# pylint: disable=missing-function-docstring

class TestTERCurveApplyTo:
    def test_adds_bg_to_source_if_source_has_no_bg(self):

        src = so._empty_sky()
        eff = eo._filter_surface()

        src = eff.apply_to(src)

        assert src.fields[-1].header["BG_SRC"] is True
        assert src.fields[-1].header["SPEC_REF"] == 1

        if PLOTS:
            for spec in src.spectra:
                wave = np.linspace(0.3, 20, 1000) * u.um
                flux = spec(wave)
                plt.semilogy(wave, flux, "b")

            eff2 = eo._filter_surface(temperature=-272)
            src = eff2.apply_to(src)

            for spec in src.spectra:
                wave = np.linspace(0.3, 20, 1000) * u.um
                flux = spec(wave)
                plt.semilogy(wave, flux, "r")
            plt.show()


class TestTERCurvePlot:
    def test_plots_only_transmission(self):
        filt = eo._filter_surface(wave_min=0.8, wave_max=2.5)
        if PLOTS:
            filt.plot(which="ter")
            plt.show()


class TestTopHatFilter:
    def test_throws_error_for_wrong_keywords(self):
        with pytest.raises(ValueError):
            tc.TopHatFilterCurve()

    def test_has_correct_profile(self):
        filt = tc.TopHatFilterCurve(transmission=0.9, wing_transmission=0.01,
                                    blue_cutoff=1., red_cutoff=2.)
        assert filt.throughput(1.*u.um).value == 0.9
        assert filt.throughput(0.998*u.um).value == 0.01


class TestDownloadableFilterCurveInit:
    def test_throws_error_if_no_keys_passed(self):
        with pytest.raises(ValueError):
            tc.DownloadableFilterCurve()

    @pytest.mark.parametrize("name_format, name",
                             [("Paranal/HAWKI.{}", "Ks"),
                              ("HST/WFC3_IR.{}", "F160W"),
                              ("JWST/NIRCam.{}", "F164N")
                              ])
    def test_returns_filter_curve_for_correct_keys(self, name, name_format):
        filt = tc.DownloadableFilterCurve(filter_name=name,
                                          filename_format=name_format)
        assert isinstance(filt, tc.DownloadableFilterCurve)
        assert isinstance(filt, tc.FilterCurve)


class TestSpanishVOFilterCurveInit:
    @pytest.mark.parametrize("observatory, instrument, filt_name",
                             [("Paranal", "HAWKI", "Ks"),
                              ("HST", "WFC3_IR", "F160W"),
                              ("JWST", "NIRCam", "F164N")])
    def test_returns_filter_as_wanted(self, observatory, instrument, filt_name):
        filt = tc.SpanishVOFilterCurve(observatory=observatory,
                                       instrument=instrument,
                                       filter_name=filt_name)
        assert isinstance(filt, tc.FilterCurve)

    def test_returns_unity_transmission_for_wrong_name(self):
        filt = tc.SpanishVOFilterCurve(observatory=None,
                                       instrument=None,
                                       filter_name=None,
                                       error_on_wrong_name=False)
        assert isinstance(filt, tc.FilterCurve)
        assert np.all([t == 1 for t in filt.data["transmission"]])


@pytest.fixture(name="fwheel", scope="class")
def _filter_wheel():
    '''Instantiate a FilterWheel'''
    return tc.FilterWheel(**{"filter_names": ["Ks", "Br-gamma"],
                             "filename_format": "TC_filter_{}.dat",
                             "current_filter": "Br-gamma"})

class TestFilterWheelInit:
    def test_initialises_correctly(self, fwheel):
        assert isinstance(fwheel, tc.FilterWheel)

    def test_loads_all_filters(self, fwheel):
        assert len(fwheel.filters) == 2

    def test_current_filter_is_filter(self, fwheel):
        assert isinstance(fwheel.current_filter, tc.FilterCurve)

    def test_current_filter_has_fov_grid_method(self, fwheel):
        assert hasattr(fwheel.current_filter, "fov_grid")

    def test_change_to_known_filter(self, fwheel):
        fwheel.change_filter('Ks')
        assert fwheel.current_filter.meta["name"] == "Ks"

    def test_change_to_unknown_filter(self, fwheel):
        with pytest.raises(ValueError):
            fwheel.change_filter('X')

    def test_plots_all_filters(self, fwheel):
        if PLOTS:
            fwheel.plot()
            plt.show()


class TestSpanishVOFilterWheelInit:
    def test_throws_exception_on_empty_input(self):
        with pytest.raises(ValueError):
            tc.SpanishVOFilterWheel()

    @pytest.mark.parametrize("observatory, instrument, default_filter",
                             [("GTC", "OSIRIS", "sdss_r_filter"),
                              ("JWST", "MIRI", "F2300C")])
    def test_returns_filter_as_wanted(self, observatory, instrument,
                                      default_filter):
        filt_wheel = tc.SpanishVOFilterWheel(observatory=observatory,
                                             instrument=instrument,
                                             current_filter=default_filter,
                                             name="test_svo_wheel")

        assert isinstance(filt_wheel, tc.FilterWheel)
        assert default_filter in filt_wheel.filters

    def test_returns_filters_with_include_str(self):
        filt_wheel = tc.SpanishVOFilterWheel(observatory="GTC",
                                             instrument="OSIRIS",
                                             current_filter="sdss_r_filter",
                                             name="test_svo_wheel",
                                             include_str="_filter")

        assert np.all(["_filter" in name for name in filt_wheel.filters])

    def test_returns_filters_with_exclude_str(self):
        filt_wheel = tc.SpanishVOFilterWheel(observatory="GTC",
                                             instrument="OSIRIS",
                                             current_filter="sdss_r_filter",
                                             name="test_svo_wheel",
                                             exclude_str="_filter")

        assert np.all(["_filter" not in name for name in filt_wheel.filters])


class TestTopHatFilterList:
    def test_throws_exception_on_empty_input(self):
        with pytest.raises(ValueError):
            tc.TopHatFilterWheel()

    def test_initialises_with_correct_input(self):
        filt_wheel = tc.TopHatFilterWheel(name="test_tophat_filter_wheel",
                                          current_filter="K",
                                          filter_names=["J", "H", "K"],
                                          transmissions=[0.9, 0.95, 0.85],
                                          wing_transmissions=[0., 0., 0.001],
                                          blue_cutoffs=[1.15, 1.45, 1.9],
                                          red_cutoffs=[1.35, 1.8, 2.4])

        assert isinstance(filt_wheel, tc.TopHatFilterWheel)
        assert filt_wheel.filters["J"].throughput(1.15*u.um) == 0.9
        assert filt_wheel.filters["J"].throughput(1.13*u.um) == 0.
        assert filt_wheel.meta["current_filter"] == "K"
