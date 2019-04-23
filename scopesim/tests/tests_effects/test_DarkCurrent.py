import numpy as np

from scopesim.detector import Detector
from scopesim.effects.electronic import DarkCurrent
from scopesim.optics.image_plane_utils import header_from_list_of_xy


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(DarkCurrent(), DarkCurrent)

    def test_initialises_with_dark_current_value_and_obs_dit_keys(self):
        dark_eff = DarkCurrent(value=0.1, OBS_DIT=10)
        assert isinstance(dark_eff, DarkCurrent)


class TestApplyTo:
    def test_applies_dark_current_with_level_of_dit(self):
        level, dit, hw = 0.5, 10, 16
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        dtcr = Detector(hdr)
        dark_eff = DarkCurrent(value=level, OBS_DIT=dit)

        dtcr = dark_eff.apply_to(dtcr)

        assert np.average(dtcr.data) == 5.0
        assert np.sum(dtcr.data) == level * dit * (hw * 2)**2




