"""Tests for METIS WCU classes"""

import pytest

from scopesim.effects.metis_wcu import BlackBodySource

@pytest.fixture(name="bbsource", scope="function")
def fixture_bbsource():
    return BlackBodySource()

class TestBlackBodySource:
    def test_initialises_correctly(self, bbsource):
        assert isinstance(bbsource, BlackBodySource)
