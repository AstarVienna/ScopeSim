# -*- coding: utf-8 -*-
"""Tests for the METIS IFU_Simple mode"""

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

import pytest
import numpy as np
from scopesim.effects.metis_ifu_simple import LineSpreadFunction

class TestLineSpreadFunction:
    def test_initialises_correctly(self):
        assert isinstance(LineSpreadFunction(lsfwidth=2.7), LineSpreadFunction)

    @pytest.mark.parametrize("width", [1.5, 2.7, 3.0, 5.2])
    def test_kernel_is_normalised(self, width):
        lsf = LineSpreadFunction(lsfwidth=width)
        assert np.isclose(lsf.kernel.sum(), 1, 0.02)


    def test_kernel_is_normalised_with_params(self):
        lsf = LineSpreadFunction(wavelen=3.0,
                                 fit_slope=3.795e-6,
                                 fit_intercept=-4.659e-7,
                                 slice_width=0.0207,
                                 pixel_scale=0.0082,
                                 spec_binwidth=1e-5)
        assert isinstance(lsf, LineSpreadFunction)
        assert np.isclose(lsf.kernel.sum(), 1, 1e-4)
