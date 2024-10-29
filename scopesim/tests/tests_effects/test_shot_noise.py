import pytest
import numpy as np

from scopesim.detector import Detector
from scopesim.effects.electronic import ShotNoise
from scopesim.optics.image_plane_utils import header_from_list_of_xy


class TestInit:
    def test_fine_when_no_keywords_are_passed(self):
        sn = ShotNoise()
        assert isinstance(sn, ShotNoise)


class TestApplyTo:
    def test_applies_dark_current_with_level_of_dit(self):
        hw = 1
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        dtcr = Detector(hdr)
        dtcr._hdu.data[0][0] = 100
        dtcr._hdu.data[0][1] = -1
        dtcr._hdu.data[1][0] = 2e33
        dtcr._hdu.data[0][1] = np.nan

        sn = ShotNoise()
        dtcr = sn.apply_to(dtcr)

        assert 1 < dtcr._hdu.data[0][0] < 1000
        assert np.isnan(dtcr._hdu.data[0][1]) or dtcr._hdu.data[0][1] == 0
        assert dtcr._hdu.data[1][0] == 2e33
        assert np.isnan(dtcr._hdu.data[0][1])
