import pytest
import numpy as np
from astropy.io.fits import Header, ImageHDU
from astropy.wcs import WCS

from scopesim.base_classes import PoorMansHeader


@pytest.fixture(scope="function")
def basic_pmh():
    hdr = ImageHDU(data=np.zeros((10, 10))).header
    return PoorMansHeader(hdr)


class TestPoorMansHeaderInit:
    def test_initialised_with_nothing(self):
        assert isinstance(PoorMansHeader(), PoorMansHeader)

    def test_initialised_with_dict(self):
        assert isinstance(PoorMansHeader({"hello": "world!"}), PoorMansHeader)

    def test_initialised_with_fits_header(self):
        hdr = ImageHDU(data=np.zeros((10, 10))).header
        assert isinstance(PoorMansHeader(hdr), PoorMansHeader)


class TestPrint:
    # FIXME: this should capture output to see if something is actually printed
    def test_print(self, basic_pmh):
        print(basic_pmh)


class TestUpdate:
    def test_updates_with_dict(self, basic_pmh):
        basic_pmh.update({"CDELT1": 0.5, "CDELT2": 1})
        assert basic_pmh["CDELT1"] == 0.5

    def test_updates_with_fits_header(self, basic_pmh):
        hdr = ImageHDU(data=np.zeros((20, 15))).header
        basic_pmh.update(hdr)
        assert basic_pmh["NAXIS1"] == 15
        assert basic_pmh["NAXIS2"] == 20

    def test_update_with_other_poormansheader(self):
        pmh1 = PoorMansHeader(ImageHDU(data=np.zeros((10, 10))).header)
        pmh2 = PoorMansHeader(ImageHDU(data=np.zeros((20, 15))).header)
        pmh1.update(pmh2)
        assert pmh1["NAXIS1"] == 15


class TestSetItem:
    def test_sets_single_value(self, basic_pmh):
        basic_pmh["CDELT1"] = 0.5
        assert basic_pmh["CDELT1"] == 0.5

    def test_sets_commented_value(self, basic_pmh):
        basic_pmh["CDELT1"] = (9001, "It's over 9000!")
        assert basic_pmh["CDELT1"] > 9000
        assert basic_pmh.comments["CDELT1"] == "It's over 9000!"


class TestIterate:
    def test_can_be_used_in_a_for_loop(self, basic_pmh):
        for key in basic_pmh:
            assert key in basic_pmh


class TestAsHeader:
    def test_returns_header_object(self, basic_pmh):
        hdr = basic_pmh.as_header()
        assert isinstance(hdr, Header)
        assert hdr["XTENSION"] == 'IMAGE'
        assert hdr.comments["XTENSION"] == 'Image extension'


class TestDict:
    def test_returns_a_dict_for_built_in_python_function_dict(self, basic_pmh):
        assert isinstance(dict(basic_pmh), dict)


class TestCompatibilityWithAstropy:
    def test_fits_header_can_update_with_it(self, basic_pmh):
        hdr = Header()
        hdr.update(basic_pmh)
        assert hdr["XTENSION"] == 'IMAGE'
        assert hdr.comments["XTENSION"] == 'Image extension'

    def test_can_make_a_wcs_from_it(self, basic_pmh):
        dic = {"CDELT1": 1, "CDELT2": 1, "CRVAL1": 5, "CRVAL2": 5,
               "CRPIX1": 0, "CRPIX2": 0, "CTYPE1": "LINEAR", "CTYPE2": "LINEAR",
               "CUNIT1": "mm", "CUNIT2": "mm", }
        basic_pmh.update(dic)
        wcs = WCS(basic_pmh)
        xw, yw = wcs.all_pix2world([10], [10], 1)
        assert isinstance(wcs, WCS)
        assert xw == 15 and yw == 15
