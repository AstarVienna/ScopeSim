# -*- coding: utf-8 -*-
"""Tests for the METIS IFU_Simple mode"""

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from scopesim.effects.metis_ifu_simple import LineSpreadFunction

class TestLineSpreadFunction:
    def test_initialises_correctly(self):
        assert isinstance(LineSpreadFunction(), LineSpreadFunction)
