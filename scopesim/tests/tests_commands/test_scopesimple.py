# -*- coding: utf-8 -*-

import pytest

from scopesim import Simulation


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestScopeSimple:
    def test_init(self):
        simple = Simulation("basic_instrument")
        assert isinstance(simple, Simulation)

    def test_instrument_name_in_str(self):
        simple = Simulation("basic_instrument")
        assert "basic_instrument" in str(simple)
