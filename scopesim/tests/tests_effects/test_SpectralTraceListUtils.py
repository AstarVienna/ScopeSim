"""Unit tests for spectral_trace_list_utils.py"""
import pytest

import numpy as np

from scopesim.effects.spectral_trace_list_utils import *

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


class TestTransform2D:
    def test_initialises_with_matrix(self):
        matrix = np.array([[1, 1], [1, 1]])
        assert isinstance(Transform2D(matrix), Transform2D)

    def test_call_gives_correct_result(self):
        x = np.random.randn()
        y = np.random.randn()

        # analytic expression
        z_a = (1 + 2 * x + 3 * x**2
               + y * (-1 - 0.5 * x + 0.3 * x**2)
               + y**2 * (0.5 - 0.3 * x + x**2))

        # matrix / Transform2D expression
        matrix = np.array([[1, 2, 3], [-1, -0.5, 0.3], [0.5, -0.3, 1]])
        td2d = Transform2D(matrix)

        assert td2d(x, y) == z_a

# ..todo: test grid=True and grid=False
