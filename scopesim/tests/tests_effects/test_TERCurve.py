import os
import pytest
from pytest import approx

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


class TestDownloadableFilterCurveInit:
    def test_throws_error_if_no_keys_passed(self):
        with pytest.raises(ValueError):
            tc.DownloadableFilterCurve()

    @pytest.mark.parametrize("name_format, name",
                             [("Paranal/HAWKI.{}", "Ks"),
                              ("HST/WFC3_IR.{}", "F160W"),
                              ("JWST/NIRCam.{}", "F164N")
                              ])
    def test_returns_fulter_curve_for_correct_keys(self, name, name_format):
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


class TestTERCurveApplyTo:
    def test_adds_bg_to_source_if_source_has_no_bg(self):

        src = so._empty_sky()
        eff = eo._filter_surface()
        src = eff.apply_to(src)

        assert src.fields[-1].header["BG_SRC"] is True
        assert src.fields[-1].header["SPEC_REF"] == 1

        wave = src.spectra[1].model.points[0]
        flux = src.spectra[1].model.lookup_table
        plt.plot(wave, flux)
        plt.show()
