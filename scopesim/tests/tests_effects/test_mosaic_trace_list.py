"""Unit tests for mosaic_trace_list.py"""

# pylint: disable=missing-function-docstring
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
import pytest

import numpy as np

from astropy.io import fits

from scopesim.utils import power_vector
from scopesim.effects.mosaic_trace_list import Transform1D

@pytest.fixture(name="tf1d", scope="class")
def fixture_tf1d():
    """Instantiate a Transform1D"""
    coeffs = np.array([2, -1, 1])
    return Transform1D(coeffs)

@pytest.fixture(name="quadratic", scope="class")
def fixture_quadratic():
    """Quadratic model, analytic and coeffients"""
    coeffs = np.array([1, -1, 2])

    def quadfunc(x):
        z_a = 1 - 1 * x + 2 * x**2
        return z_a

    def dquad_dx(x):
        return -1 + 4 * x

    return {'coeffs': coeffs,
            'function': quadfunc,
            'gradient': dquad_dx}

class TestTransform1D:
    """Tests for Transform1D()"""
    def test_initialises_with_coeffs(self, tf1d):
        assert isinstance(tf1d, Transform1D)

    def test_call_gives_correct_result(self, quadratic):
        x = np.random.randn()

        # coefficients and explicit function
        tf1d = Transform1D(quadratic['coeffs'])
        assert tf1d(x) == quadratic['function'](x)

    def test_gradient_gives_correct_result(self, quadratic):
        x = np.random.randn()

        tf2d = Transform1D(quadratic['coeffs'])
        tf2d_grad = tf2d.gradient()

        assert tf2d_grad(x) == quadratic['gradient'](x)

    def test_fit_gives_correct_coeffs(self):
        x = np.linspace(0, 1, 10)
        y = 1 - 0.5 * x + 2.3 * x**2 - 3 * x**3

        coeffs = np.array([1, -0.5, 2.3, -3])
        tf1d = Transform1D.fit(x, y, degree=3)

        assert tf1d.coeffs == pytest.approx(coeffs)
