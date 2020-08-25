import os

import numpy as np
import pytest
from astropy.table import Table

from scopesim import rc
from scopesim.effects import Effect, SurfaceList

from scopesim.tests.mocks.py_objects import effects_objects as eo


MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]


@pytest.fixture()
def surf_list_file():
    fname = os.path.join(MOCK_PATH, "LIST_mirrors_MICADO_Wide.tbl")
    return fname


class TestEffectInit:
    def test_initialises_with_no_input(self):
        assert isinstance(Effect(), Effect)

    def test_initalising_with_arrays_creates_table(self):
        eff = Effect(array_dict={"x": [-1, 0, 1], "y": [1, 0, 1],
                                 "flux": [1, 2, 3]})
        assert isinstance(eff, Effect)
        assert np.sum(eff.table["flux"]) == 6

    def test_has_method_apply_to(self):
        assert hasattr(Effect(), "apply_to")

    def test_has_method_waveset(self):
        assert hasattr(Effect(), "fov_grid")


class TestEffectReport:
    def test_report_returns_full_rst_text(self):
        det_list = eo._detector_list()
        det_list.meta.update(
            {"report_plot_include": False,
             "report_table_caption": "The dimensions of the MICADO central detector"})
        rst_str = det_list.report()
        assert "MICADO H4RG-15 FPA" in rst_str
        assert "E-MCD-FPA-572089EB.uda" in rst_str
        assert "The dimensions of the MICADO central detector" in rst_str


@pytest.mark.usefixtures("surf_list_file")
class TestSurfaceListInit:
    def test_initialises_with_nothing(self):
        assert isinstance(SurfaceList(), SurfaceList)

    def test_initialises_with_valid_filename(self, surf_list_file):
        surf_list = SurfaceList(filename=surf_list_file)
        assert isinstance(surf_list, SurfaceList)
        assert isinstance(surf_list.data, Table)


