import pytest
import os
import numpy as np

from scopesim import rc
from scopesim.effects import LinearityCurve

from scopesim.tests.mocks.py_objects.detector_objects import _basic_detector

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))

for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


class TestInit:
    def test_initialises_with_direct_ndit(self):
        assert isinstance(LinearityCurve(ndit=1), LinearityCurve)

    # Not relevant because the bangkeys are now called just-in-time
    # def test_initialises_with_bangkey_ndit(self):
    #     rc.__currsys__["!OBS.ndit"] = 2
    #     lincurve = LinearityCurve(ndit="!OBS.ndit")
    #     assert lincurve.meta["ndit"] == 2

    def test_throws_error_with_no_keywords(self):
        with pytest.raises(ValueError):
            LinearityCurve()

    def test_initialises_with_filename(self):
        lincurve = LinearityCurve(ndit=2, filename="test_linearity.dat")
        assert "incident" in lincurve.table.colnames
        assert "measured" in lincurve.table.colnames

    def test_initialises_with_arrays(self):
        lincurve = LinearityCurve(array_dict={"incident": [0, 50, 100],
                                              "measured": [0, 75, 100]},
                                  ndit=2)
        assert "incident" in lincurve.table.colnames
        assert "measured" in lincurve.table.colnames

    def test_initialises_with_vectors(self):
        lincurve = LinearityCurve(incident=[0, 50, 100],
                                  measured=[0, 75, 100],
                                  ndit=2)
        assert "incident" in lincurve.meta
        assert "measured" in lincurve.meta

class TestApplyTo:
    @pytest.mark.parametrize("in_flux, out_flux", [(20, 10), (45, 52.5),
                                                   (60, 60), (1e5, 60)])
    def test_only_applied_to_detector(self, in_flux, out_flux):
        dtcr = _basic_detector()
        dtcr._hdu.data[0, 0] = in_flux

        lincurve = LinearityCurve(ndit=1, filename="test_linearity.dat")
        new_dtcr = lincurve.apply_to(dtcr)
        assert new_dtcr._hdu.data[0, 0] == out_flux
        assert np.min(new_dtcr._hdu.data) == 0

    def test_bypasses_non_detector_objects(self):
        lincurve = LinearityCurve(ndit=1, filename="test_linearity.dat")
        new_dict = lincurve.apply_to({"gigawatts": 1.21})
        assert new_dict["gigawatts"] == 1.21
