"""Unit tests for module scopesim.utils"""

import os
import logging
import pytest
import numpy as np
from astropy import wcs
from astropy.io import ascii as ioascii, fits
from astropy.table import Table

from scopesim import utils

from scopesim import rc


class TestFindFile:
    """Tests of function scopesim.utils.find_file"""

    def test_fails_if_filename_not_a_string(self):
        # python 3.6: TypeError
        # python 3.4, 3.5: AttributeError (change in os.path.isabs)
        with pytest.raises((TypeError, AttributeError)):
            utils.find_file(1.2, rc.__search_path__)

    def test_passes_if_file_exists(self):
        filename = 'utils.py'
        assert utils.find_file(filename, rc.__search_path__)

    @pytest.mark.parametrize("throw_error", [True, False])
    def test_throws_error_if_file_doesnt_exist(self, throw_error):
        rc.__currsys__["!SIM.file.error_on_missing_file"] = throw_error
        filename = 'utils987654.pz'

        if throw_error:
            with pytest.raises(ValueError):
                utils.find_file(filename, rc.__search_path__)
        else:
            assert utils.find_file(filename, rc.__search_path__) is None

        rc.__currsys__["!SIM.file.error_on_missing_file"] = False

    def test_ignores_none_objects_in_search_path_list(self):
        filename = 'utils.py'
        new_filename = utils.find_file(filename, [None] + rc.__search_path__)
        assert filename in new_filename


class TestAirmassZendist:
    """Tests conversion between airmass and zenith distance"""

    def test_airmass2zendist_pass_for_known_quanities_AM_1_equals_ZD_0(self):
        assert np.allclose(utils.airmass2zendist(1.0), 0)

    def test_pass_for_known_quanities_AM_2_equals_ZD_sqrt2(self):
        assert np.allclose(utils.airmass2zendist(np.sqrt(2)), 45)

    def test_zendist2airmass_pass_for_known_quanities_ZD_0_equals_AM_1(self):
        assert np.allclose(utils.zendist2airmass(0), 1.0)

    def test_zendist2airmass_pass_for_known_quanities_ZD_60_equals_ZD_2(self):
        assert np.allclose(utils.zendist2airmass(60), 2.0)

    def test_zendist2airmass_undoes_exactly_what_airmass2zendist_does(self):
        airmass = 1.78974234
        assert np.allclose(utils.zendist2airmass(utils.airmass2zendist(airmass)),
                           airmass)

    def test_airmass2zendist_undoes_exactly_what_zendist2airmass_does(self):
        zendist = 12.31334
        assert np.allclose(utils.airmass2zendist(utils.zendist2airmass(zendist)),
                           zendist)


class TestParallacticAngle:
    """Tests of function scopesim.utils.parallactic_angle"""

    def test_parallactic_angle_negative_east_of_meridian(self):
        assert utils.parallactic_angle(-1, 0, -24) < 0

    def test_parallactic_angle_positive_west_of_meridian(self):
        assert utils.parallactic_angle(1, 0, -24) > 0

    def test_parallactic_angle_zero_on_meridian(self):
        assert utils.parallactic_angle(0, 0, 24) == 0

    def test_specific_example_from_Ball_1908(self):
        """Test: Example from Ball (1908), p.92"""
        ha = -3.                 # 3 hours east
        de = 38 + 9/60.          # decl 38d09m
        lat = 53 + 23/60.        # lat  53d23m
        eta0 = - (48 + 41/60.)   # result -48d41m

        eta = utils.parallactic_angle(ha, de, lat)

        # should agree to within 1 arcmin
        assert np.allclose(eta, eta0, atol=1/60.)

    def test_setting_object_on_the_equator_is_90_minus_latitude(self):
        """
        For a setting object on the equator, the parallactic angle
        is 90 - lat
        """
        lat = np.random.rand(10) * 180 - 90
        pa = utils.parallactic_angle(6, 0, lat)

        assert np.allclose(pa, 90. - lat)


class TestDerivPolynomial2D:
    """Tests of scopesim.utils.deriv_polynomial2d"""

    def test_derivative_of_2D_polynomial_equal_to_analytical_derivative(self):
        from astropy.modeling.models import Polynomial2D

        ximg, yimg = np.meshgrid(np.linspace(-1, 1, 101),
                                 np.linspace(-1, 1, 101))
        poly = Polynomial2D(2, c0_0=1, c1_0=2, c2_0=3,
                            c0_1=-1.5, c0_2=0.4, c1_1=-2)
        # Expected values
        y_x = 2 + 6 * ximg - 2 * yimg
        y_y = -1.5 + 0.8 * yimg - 2 * ximg

        dpoly_x, dpoly_y = utils.deriv_polynomial2d(poly)
        # Computed values
        y_x_test = dpoly_x(ximg, yimg)
        y_y_test = dpoly_y(ximg, yimg)

        assert np.allclose(y_x, y_x_test)
        assert np.allclose(y_y, y_y_test)


class TestConvertCommentsToDict:
    def test_converts_list_of_strings_to_dict_if_comments_in_table_meta(self):
        tbl = ioascii.read("""# key1 : val 1 
                              # key2 : extra long entry
                              col1    col2
                              0       1 """)
        dic = utils.convert_table_comments_to_dict(tbl)
        assert dic["key1"] == "val 1"
        assert len(dic) == 2

    def test_returns_empty_dict_if_comments_not_in_table_meta(self):
        tbl = ioascii.read("""col1    col2
                              0       1 """)
        dic = utils.convert_table_comments_to_dict(tbl)
        assert dic == {}

    def test_returns_input_if_conversion_doesnt_work(self):
        tbl_str = """
        # key1 : val 1 
        # 
        # key2
        col1    col2
        0       1 """
        tbl = ioascii.read(tbl_str)
        dic = utils.convert_table_comments_to_dict(tbl)
        assert dic == tbl.meta["comments"]


class TestHasWcsKeys:
    def test_fails_if_header_does_not_have_all_keys(self):
        assert not utils.has_needed_keywords(fits.Header())

    def test_passes_if_header_does_have_all_keys(self):
        hdr = wcs.WCS().to_header()
        hdr["NAXIS1"] = 100
        assert utils.has_needed_keywords(hdr)

    def test_passes_if_header_does_have_all_keys_and_suffix(self):
        hdr = wcs.WCS(key="D").to_header()
        hdr["NAXIS1"] = 100
        assert utils.has_needed_keywords(hdr, "D")


class TestFromCurrSys:
    def test_converts_string(self):
        assert utils.from_currsys("!SIM.random.seed") is None

    def test_converts_list(self):
        assert utils.from_currsys(["!SIM.random.seed"]*3)[2] is None

    def test_converts_numpy_array(self):
        assert utils.from_currsys(np.array(["!SIM.random.seed"]*2))[1] is None

    def test_converts_dict(self):
        assert utils.from_currsys({"seed": "!SIM.random.seed"})["seed"] is None

    def test_converts_astropy_table(self):
        tbl = Table(data=[["!SIM.random.seed"]*2, ["!SIM.random.seed"]*2],
                    names=["seeds", "seeds2"])
        assert utils.from_currsys(tbl["seeds2"][1]) is None


class TestSetupLoggers:
    def test_no_loggers_present_in_initial_scopesim_load(self):
        hdlrs = logging.getLogger().handlers

        assert not any([hdlr.name == "scopesim_file_logger" for hdlr in hdlrs])

    def test_setup_loggers_has_no_loggers(self):
        utils.setup_loggers(log_to_file=False, log_to_console=False)
        hdlrs = logging.getLogger().handlers

        assert not any([hdlr.name == "scopesim_file_logger" for hdlr in hdlrs])

    def test_log_file_exist_when_file_logging_turned_on(self):
        utils.setup_loggers(log_to_file=True)
        filepath = rc.__config__["!SIM.logging.file_path"]
        level = rc.__config__["!SIM.logging.file_level"]
        f_handler = logging.getLogger().handlers[-1]

        assert os.path.exists(filepath)
        assert f_handler.level == logging._checkLevel(level)

        logging.shutdown()
        os.remove(filepath)

    def test_console_logger_has_output(self, capsys):
        utils.setup_loggers(log_to_console=True)
        utils.set_logger_level("console", "INFO")
        logging.info("Hello World!")
        captured = capsys.readouterr()
        assert "Hello World!" in captured.out

        logging.shutdown()
