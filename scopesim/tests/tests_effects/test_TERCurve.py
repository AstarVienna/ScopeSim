import pytest
from pytest import approx

import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u

from scopesim.effects import ter_curves as tc


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

