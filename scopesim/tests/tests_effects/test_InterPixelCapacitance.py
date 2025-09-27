# -*- coding: utf-8 -*-
"""Tests for InterPixelCapacitance effect."""

import pytest
import numpy as np

from scopesim.effects.electronic import InterPixelCapacitance as IPC

# pylint: disable=missing-class-docstring,missing-function-docstring
class TestInit:
    def test_initialises_correctly(self):
        ipc = IPC()
        assert isinstance(ipc, IPC)

    def test_initialised_with_nothing(self):
        ipc = IPC()
        assert np.all(ipc.kernel == np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]]))

    def test_initialises_with_kernel(self):
        kern = np.random.rand(3, 3)
        ipc = IPC(kernel=kern)
        assert np.all(ipc.kernel == kern)

    @pytest.mark.parametrize(
        "a_edge, a_corner, a_cross, kern",
        [(0, 0, 0, [[0, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]]),
         (0.02, 0, 0, [[0, 0.02, 0],
                       [0.02, 0.92, 0.02],
                       [0, 0.02, 0]]),
         (0.02, 0.002, 0, [[0.002, 0.02, 0.002],
                           [0.02, 0.912, 0.02],
                           [0.002, 0.02, 0.002]]),
         (0.0145, 0.0011, 0.0018, [[0.0011, 0.0127, 0.0011],
                                   [0.0163, 0.9376, 0.0163],
                                   [0.0011, 0.0127, 0.0011]])
         ])
    def test_initialises_with_params(self, a_edge, a_corner, a_cross, kern):
        ipc = IPC(alpha_edge=a_edge, alpha_corner=a_corner, alpha_cross=a_cross)
        assert np.allclose(ipc.kernel, np.asarray(kern))
