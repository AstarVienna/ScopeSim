# -*- coding: utf-8 -*-
"""To run these tests in a script, set:
sim.rc.__config__["!SIM.file.local_packages_path"] = "./scopesim/tests/mocks/"
"""

import pytest

from scopesim import Simulation
from scopesim.source import source_templates as st


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestScopeSimple:
    def test_init(self):
        simple = Simulation("basic_instrument")
        assert isinstance(simple, Simulation)

    def test_instrument_name_in_str(self):
        simple = Simulation("basic_instrument")
        assert "basic_instrument" in str(simple)

    def test_default_mode_works(self):
        simple = Simulation("basic_instrument")
        assert simple.mode == "imaging"

    def test_init_mode_works(self):
        simple = Simulation("basic_instrument", "spectroscopy")
        assert simple.mode == "spectroscopy"

    def test_full_workflow_runs(self):
        simple = Simulation("basic_instrument")
        src = st.star(flux=15)
        out1 = simple(src, dit=1, ndit=1)
        out2 = simple.readout(dit=1e5, ndit=1)
        assert out2[1].data.max() > out1[1].data.max()

    def test_plotting_runs(self):
        simple = Simulation("basic_instrument")
        src = st.star(flux=9)
        simple(src, dit=10, ndit=1)
        fig, ax = simple.plot(adjust_scale=True)
        assert ax is not None
