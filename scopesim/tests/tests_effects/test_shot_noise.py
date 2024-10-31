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
        """Test whether ShotNoise is applied correctly."""
        hw = 1
        # A normal value that works with numpy.random.poisson.
        value_normal = 1000
        # A value too large fro numpy.random.poisson.
        value_too_high = 1e20
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        dtcr = Detector(hdr)
        dtcr._hdu.data[0][0] = value_normal
        dtcr._hdu.data[0][1] = -1
        dtcr._hdu.data[1][0] = value_too_high
        dtcr._hdu.data[1][1] = np.nan

        sn = ShotNoise()
        dtcr = sn.apply_to(dtcr)

        # Ensure that the values have changed.
        assert dtcr._hdu.data[0][0] != value_normal
        assert dtcr._hdu.data[1][0] != value_too_high

        # Sensibility checks on the values.
        assert 1 < dtcr._hdu.data[0][0] < value_normal * 2
        assert np.isnan(dtcr._hdu.data[0][1]) or dtcr._hdu.data[0][1] == 0
        assert np.isnan(dtcr._hdu.data[1][1])
