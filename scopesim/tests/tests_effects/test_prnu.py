from matplotlib import pyplot as plt
import numpy as np
import pytest

from scopesim.detector import Detector
from scopesim.effects.electronic import PixelResponseNonUniformity
from scopesim.optics.image_plane_utils import header_from_list_of_xy

PLOTS = False

def make_detector(value=1000, size=10):
    hdr = header_from_list_of_xy([-size/2, size/2], [-size/2, size/2], 1, "D")
    dtcr = Detector(hdr)
    dtcr._hdu.data[:] = value
    return dtcr


class TestApplyTo:
    def test_output_std_matches_prnu_std(self):
        """Relative std of output should be close to prnu_std."""
        prnu_std = 0.05
        dtcr = make_detector(value=1000, size=100)
        PixelResponseNonUniformity(prnu_std=prnu_std, prnu_seed=42).apply_to(dtcr)
        rel_std = dtcr._hdu.data.std() / dtcr._hdu.data.mean()
        assert abs(rel_std - prnu_std) < 0.01

    def test_gain_map_is_reused(self):
        """Same map should be applied on repeated calls (fixed pattern)."""
        dtcr1 = make_detector()
        dtcr2 = make_detector()
        prnu = PixelResponseNonUniformity(prnu_std=0.01, prnu_seed=42)
        prnu.apply_to(dtcr1)
        prnu.apply_to(dtcr2)
        np.testing.assert_array_equal(dtcr1._hdu.data, dtcr2._hdu.data)

    def test_dict_prnu_std(self):
        """Dict mode: different amplitude per detector ID."""
        hdr = header_from_list_of_xy([-5, 5], [-5, 5], 1, "D")
        dtcr = Detector(hdr)
        dtcr.meta["id"] = "H2RG"
        dtcr._hdu.data[:] = 1000
        prnu = PixelResponseNonUniformity(
            prnu_std={"H2RG": 0.005, "GeoSnap": 0.020}, prnu_seed=42)
        prnu.apply_to(dtcr)
        assert dtcr._hdu.data.std() > 0

    def test_multiplicative_zero_signal(self):
        """Zero signal should remain zero — confirms multiplicative (not additive)."""
        dtcr = make_detector(value=0)
        PixelResponseNonUniformity(prnu_std=0.01, prnu_seed=42).apply_to(dtcr)
        assert dtcr._hdu.data.sum() == 0

    def test_non_detector_passthrough(self):
        """Non-Detector objects should be returned unchanged."""
        prnu = PixelResponseNonUniformity(prnu_std=0.01)
        result = prnu.apply_to("not a detector")
        assert result == "not a detector"

    def test_invalid_prnu_std_raises(self):
        prnu = PixelResponseNonUniformity(prnu_std="not a number nor a dict")
        with pytest.raises(TypeError):
            prnu.apply_to(make_detector())

    def test_plot_raises_before_simulation(self):
        prnu = PixelResponseNonUniformity(prnu_std=0.01)
        with pytest.raises(RuntimeError):
            prnu.plot()

    def test_plot_returns_figure(self):
        prnu = PixelResponseNonUniformity(prnu_std=0.01, prnu_seed=42)
        prnu.apply_to(make_detector())
        if PLOTS:
            assert isinstance(prnu.plot(), plt.Figure)
