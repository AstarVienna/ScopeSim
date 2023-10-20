import pytest
import numpy as np

from scopesim.effects import LinearityCurve

from scopesim.tests.mocks.py_objects.detector_objects import _basic_detector


@pytest.fixture
def basic_lincurve(mock_path):
    fname = str(mock_path / "test_linearity.dat")
    return LinearityCurve(ndit=1, filename=fname)


class TestInit:
    def test_initialises_with_direct_ndit(self):
        assert isinstance(LinearityCurve(ndit=1), LinearityCurve)

    def test_throws_error_with_no_keywords(self):
        with pytest.raises(ValueError):
            LinearityCurve()

    def test_initialises_with_filename(self, mock_path):
        lincurve = LinearityCurve(
            ndit=2, filename=str(mock_path / "test_linearity.dat"))
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
    def test_only_applied_to_detector(self, in_flux, out_flux, basic_lincurve):
        dtcr = _basic_detector()
        dtcr._hdu.data[0, 0] = in_flux

        new_dtcr = basic_lincurve.apply_to(dtcr)
        assert new_dtcr._hdu.data[0, 0] == out_flux
        assert np.min(new_dtcr._hdu.data) == 0

    def test_bypasses_non_detector_objects(self, basic_lincurve):
        new_dict = basic_lincurve.apply_to({"gigawatts": 1.21})
        assert new_dict["gigawatts"] == 1.21
