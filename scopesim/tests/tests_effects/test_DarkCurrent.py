import pytest
import numpy as np

from scopesim.detector import Detector
from scopesim.effects.electronic import DarkCurrent
from scopesim.optics.image_plane_utils import header_from_list_of_xy


class TestInit:
    def test_fails_when_no_keywords_are_passed(self):
        with pytest.raises(ValueError):
            DarkCurrent()

    def test_initialises_with_dark_current_value_and_obs_dit_keys(self):
        dark_eff = DarkCurrent(value=0.1, dit=10, ndit=1)
        assert isinstance(dark_eff, DarkCurrent)


class TestApplyTo:
    def test_applies_dark_current_with_level_of_dit(self):
        level, dit, hw = 0.5, 10, 16
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        dtcr = Detector(hdr)
        dtcr.meta["dit"] = dit
        dtcr.meta["ndit"] = 1
        dark_eff = DarkCurrent(value=level, **dtcr.meta)

        dtcr = dark_eff.apply_to(dtcr)

        assert np.average(dtcr.data) == 5.0
        assert np.sum(dtcr.data) == level * dit * (hw * 2)**2


