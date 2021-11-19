"""Unit tests for spectral_trace_list_utils.py"""

# pylint: disable=no-self-use
# pylint: disable=missing-function-docstring
# pylint: disable=invalid-name

import pytest

import numpy as np

from scopesim.effects.spectral_trace_list_utils import Transform2D, power_vector

class TestPowerVec:
    """Test function power_vector()"""
    def test_gives_correct_result(self):
        res = power_vector(-2, 4)
        assert all(a == b for a, b in zip(res, [1, -2, 4, -8, 16]))

    def test_fails_with_negative_degree(self):
        with pytest.raises(ValueError):
            power_vector(np.pi, -2)

    def test_fails_with_non_integer_degree(self):
        with pytest.raises(ValueError):
            power_vector(np.pi, np.pi)


@pytest.fixture(name="tf2d", scope="class")
def fixture_tf2d():
    """Instantiate a Transform2D"""
    matrix = np.array([[1, 1], [0, 1]])
    return Transform2D(matrix)

@pytest.fixture(name="quadratic", scope="class")
def fixture_quadratic():
    """Quadratic model, analytic and matrix
    """
    matrix = np.array([[1, 2, 3], [-1, -0.5, 0.3], [0.5, -0.3, 1]])

    def quadfunc(x, y):
        z_a = (1 + 2 * x + 3 * x**2
               + y * (-1 - 0.5 * x + 0.3 * x**2)
               + y**2 * (0.5 - 0.3 * x + x**2))
        return z_a

    def dquad_dx(x, y):
        return (2 + 6 * x + y * (-0.5 + 0.6 * x)
                + y**2 * (-0.3 + 2 * x))

    def dquad_dy(x, y):
        return (-1 - 0.5 * x + 0.3 * x**2
                + 2 * y * (0.5 - 0.3 * x + x**2))

    return {'matrix': matrix,
            'function': quadfunc,
            'gradient': (dquad_dx, dquad_dy)}

class TestTransform2D:
    """Tests for Transform2D()"""
    def test_initialises_with_matrix(self, tf2d):
        assert isinstance(tf2d, Transform2D)

    def test_call_gives_correct_result(self, quadratic):
        x = np.random.randn()
        y = np.random.randn()

        # matrix and explicit function
        td2d = Transform2D(quadratic['matrix'])
        assert td2d(x, y) == quadratic['function'](x, y)

    def test_gradient_gives_correct_result(self, quadratic):
        x = np.random.randn()
        y = np.random.randn()

        tf2d = Transform2D(quadratic['matrix'])
        tf2d_grad = tf2d.gradient()

        assert tf2d_grad[0](x, y) == quadratic['gradient'][0](x, y)
        assert tf2d_grad[1](x, y) == quadratic['gradient'][1](x, y)

    def test_grid_true_gives_correct_shape(self, tf2d):
        n_x, n_y = 12, 3
        res = tf2d(np.ones(n_x), np.ones(n_y), grid=True)
        assert res.shape == (n_y, n_x)

    def test_grid_false_give_correct_shape(self, tf2d):
        npts = 12
        res = tf2d(np.ones(npts), np.ones(npts), grid=False)
        assert res.shape == (npts,)

    def test_grid_false_fails_with_unequal_lengths(self, tf2d):
        n_x, n_y = 12, 3
        with pytest.raises(ValueError):
            tf2d(np.ones(n_x), np.ones(n_y), grid=False)

    def test_fit_gives_correct_matrix(self):
        xx, yy = np.meshgrid(np.arange(5), np.arange(5))
        zz = 1. + xx - yy

        matrix = np.array([[1, 1], [-1, 0]])
        tf2d = Transform2D.fit(xx, yy, zz, degree=1)

        assert tf2d.matrix == pytest.approx(matrix)

    def test_grid_false_shape_is_preserved(self, tf2d):
        n_x, n_y = 4, 2
        res = tf2d(np.ones((n_y, n_x)), np.ones((n_y, n_x)), grid=False)
        assert res.shape == (n_y, n_x)
