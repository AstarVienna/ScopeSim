import pytest
from pytest import raises
from astropy.io import fits
import numpy as np

from scopesim.effects import ExtraFitsKeywords
from scopesim.effects import fits_headers as fh


@pytest.fixture(scope="function")
def comb_hdul():
    pri = fits.PrimaryHDU(header=fits.Header({"EXTNAME": "PriHDU"}))
    im = fits.ImageHDU(header=fits.Header({"EXTNAME": "ImHDU"}))
    tbl = fits.BinTableHDU(header=fits.Header({"EXTNAME": "BinTblHDU"}))
    hdul = fits.HDUList([pri, im, tbl])

    return hdul


@pytest.mark.usefixtures("comb_hdul")
class TestGetRelevantExtensions:
    def test_works_for_ext_name(self, comb_hdul):
        dic = {"ext_name": "PriHDU"}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [0]

        assert np.all([ans in exts for ans in answer])
        assert len(exts) == len(answer)

    def test_works_for_ext_number(self, comb_hdul):
        dic = {"ext_number": [1, 2, 3]}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [1, 2]

        assert np.all([ans in exts for ans in answer])
        assert len(exts) == len(answer)

    def test_works_for_ext_type(self, comb_hdul):
        dic = {"ext_type": ["ImageHDU", "PrimaryHDU"]}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [0, 1]

        assert np.all([ans in exts for ans in answer])
        assert len(exts) == len(answer)


class TestFlattenDict:
    def test_works(self):
        dic = {"HIERARCH":
                   {"ESO":
                        {"ATM":
                             {"PWV": 1.0, "AIRMASS": 2.0},
                         "DPR": {"TYPE": "DARK"}},
                    "SIM": {"area": "!TEL.area"}
                    }
               }
        flat_dict = fh.flatten_dict(dic)
        assert flat_dict["HIERARCH ESO ATM PWV"] == 1.0
        assert flat_dict["HIERARCH ESO ATM AIRMASS"] == 2.0
        assert flat_dict["HIERARCH ESO DPR TYPE"] == "DARK"
        assert flat_dict["HIERARCH SIM area"] == "!TEL.area"
