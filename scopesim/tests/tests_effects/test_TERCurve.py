
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


def _filter_wheel():
    wheel = tc.FilterWheel(**{"filter_names": ["Ks", "Br-gamma"],
                              "filename_format": "TC_filter_{}.dat",
                              "current_filter": "Br-gamma"})
    return wheel


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


class TestFilterWheelInit:
    def test_loads_all_filters(self):
        wheel = _filter_wheel()
        assert len(wheel.filters) == 2

    def test_current_filter_can_be_changed(self):
        wheel = _filter_wheel()
        assert wheel.current_filter.meta["name"] == "Br-gamma"

        wheel.meta["current_filter"] = "Ks"
        assert wheel.current_filter.meta["name"] == "Ks"

    def test_plots_all_filters(self):
        wheel = _filter_wheel()
        if PLOTS:
            wheel.plot()
            plt.show()


