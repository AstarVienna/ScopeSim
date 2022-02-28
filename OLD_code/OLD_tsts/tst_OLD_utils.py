"""
Testing the functions behind a UserCommands object
"""
import pytest
import logging
import os

from scopesim import rc
from OLD_code import OLD_user_commands_utils as cmd_utils


class TestParseConfig:
    def test_empty_str_returns_empty_dict(self):
        dic = cmd_utils.lines_to_dict([])
        assert len(dic) == 0

    def test_correctly_formatted_string_returns_ordered_dict(self):
        valid_lines = ["NO_COMMENT    True", "WITH_COMMENT    foo   # bar"]
        dic = cmd_utils.lines_to_dict(valid_lines)
        assert dic["NO_COMMENT"] is True
        assert dic["WITH_COMMENT"] == "foo"

    def test_bogus_string_raises_unknown_response(self):
        invalid_lines = ["One point twenty-one Jigawatts!"]
        dic = cmd_utils.lines_to_dict(invalid_lines)

        logging.warning("Bogus data passes silently through parse_config()")
        pass

    def test_raises_exception_if_input_is_not_list(self):
        invalid_lines = "One point twenty-one Jigawatts!"
        with pytest.raises(ValueError):
            dic = cmd_utils.lines_to_dict(invalid_lines)

    def test_correctly_coverts_none_string_to_none_type(self):
        dic = cmd_utils.lines_to_dict(["A None"])
        assert dic["A"] is None

    def test_correctly_coverts_true_string_to_true_type(self):
        dic = cmd_utils.lines_to_dict(["A True"])
        assert dic["A"] is True

    def test_correctly_coverts_number_string_to_float_type(self):
        dic = cmd_utils.lines_to_dict(["A 7"])
        assert dic["A"] == 7.


class TestReadConfig:
    def test_passes_for_correctly_formatted_multiline_string(self):
        multiline_string = """
        LINE_ONE   hello  # this is the first line
        LINE_TWO  world!  # this is the second line
        """
        dic = cmd_utils.read_config(multiline_string)
        assert dic["LINE_ONE"] == "hello"
        assert dic["LINE_TWO"] == "world!"

    def test_raises_exception_for_incorrectly_formatted_multiline_string(self):
        dodgy_multiline_string = """
        LINE_ONE  # this lines is missing a value
        """
        with pytest.raises(ValueError):
            cmd_utils.read_config(dodgy_multiline_string)

    def test_passes_when_given_filename_that_exist(self):
        rc_file = os.path.join(rc.__pkg_dir__, "OLD.scopesimrc")
        rc_dict = cmd_utils.read_config(rc_file)
        assert rc_dict["SIM_SUB_PIXEL_FLAG"] is False

    def test_raises_exception_if_filename_doesnt_exist(self):
        with pytest.raises(ValueError):
            cmd_utils.read_config("bogus.txt")

    def test_raise_exception_if_input_is_not_string(self):
        with pytest.raises(ValueError):
            cmd_utils.read_config(["hello", "world"])


class TestUpdateConfig:
    def test_update_happens_when_given_correctly_formatted_string(self):
        orig_dic = {"A" : "B"}
        update_str = "A 2"
        new_dic = cmd_utils.update_config(update_str, orig_dic)
        assert new_dic["A"] == 2
        assert orig_dic["A"] == 2

    def test_update_happens_when_given_other_dict(self):
        orig_dic = {"A" : "B"}
        update_dic = {"A" : 2}
        new_dic = cmd_utils.update_config(update_dic, orig_dic)
        assert new_dic["A"] == 2
        assert orig_dic["A"] == 2

    def test_raise_exception_when_not_given_string_or_dict(self):
        orig_dic = {"A": "B"}
        with pytest.raises(ValueError):
            cmd_utils.update_config(["a", "b"], orig_dic)


class TestConvertDictStringsToPythonTypes:
    def test_nones_are_replaces(self):
        dic = {"A" : "none", "B" : "None", "C" : "NONE", "D" : "nOnE"}
        new_dic = cmd_utils.convert_dict_strings_to_python_types(dic)
        for key in dic:
            assert dic[key] is None


class TestStrToPythonType:
    def test_conversion_works(self):
        assert cmd_utils.str_to_python_type("none") is None
        assert cmd_utils.str_to_python_type("True") is True
        assert cmd_utils.str_to_python_type("42") == 42.
        assert cmd_utils.str_to_python_type("Geronimo!") == "Geronimo!"

    def test_throws_error_if_input_is_not_string(self):
        assert cmd_utils.str_to_python_type(True) is True
        assert cmd_utils.str_to_python_type(None) is None
        assert cmd_utils.str_to_python_type(42) is 42
