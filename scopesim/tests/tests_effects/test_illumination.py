"""Tests for illumination effects"""
import pytest
import numpy as np

from scopesim.effects.illumination import Illumination, gaussian2d, quadratic_vignetting
from scopesim.optics.image_plane import ImagePlane
from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square

PLOTS = False


@pytest.fixture
def imageplane():
    ip = ImagePlane(_image_hdu_square().header)
    ip.hdu.data = np.ones((100, 100), dtype=np.float64)
    return ip


# --- gaussian2d ---

def test_gaussian2d_peak_at_centre():
    result = gaussian2d((101, 101))
    assert result[50, 50] == pytest.approx(1.0)

def test_gaussian2d_values_leq_amp():
    result = gaussian2d((101, 101))
    assert result.max() <= 1.0 + 1e-12


# --- quadratic_vignetting ---

def test_quadratic_vignetting_centre_is_one():
    result = quadratic_vignetting((101, 101))
    assert result[50, 50] == pytest.approx(1.0)

def test_quadratic_vignetting_values_in_range():
    result = quadratic_vignetting((101, 101))
    assert np.all(result >= 0.0) and np.all(result <= 1.0)


# --- Illumination ---

def test_illumination_instantiates():
    assert isinstance(Illumination(), Illumination)

def test_illumination_apply_to_returns_imageplane(imageplane):
    eff = Illumination()
    assert eff.apply_to(imageplane) is imageplane

def test_illumination_apply_to_skips_non_imageplane():
    eff = Illumination()
    obj = object()
    assert eff.apply_to(obj) is obj

def test_illumination_modifies_data(imageplane):
    original = imageplane.hdu.data.copy()
    Illumination().apply_to(imageplane)
    assert not np.array_equal(imageplane.hdu.data, original)

def test_illumination_caches_map(imageplane):
    eff = Illumination()
    eff.apply_to(imageplane)
    assert eff._map is not None and eff._map_shape == (100, 100)

def test_illumination_make_map_shape_and_dtype():
    eff = Illumination()
    m = eff._make_map((80, 60))
    assert m.shape == (80, 60) and m.dtype == np.float32

def test_illumination_plot_raises_before_apply():
    with pytest.raises(RuntimeError):
        Illumination().plot()

