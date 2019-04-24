import sys
import os

from scopesim import __version__, rc
from scopesim.commands.user_commands_utils import read_config


class TestBasicLoading:

    def test_pkg_dir_in_search_path(self):
        assert rc.__pkg_dir__ in rc.__search_path__

    def test_data_dir_in_search_path(self):
        assert rc.__data_dir__ in rc.__search_path__

    def test_data_dir_in_pkg_dir(self):
        if sys.version_info >= (3,5):
            cpath = os.path.commonpath([rc.__pkg_dir__, rc.__data_dir__])
            assert cpath == rc.__pkg_dir__
        else:
            cpath = os.path.commonprefix([rc.__pkg_dir__, rc.__data_dir__])
            assert cpath == rc.__pkg_dir__

    def test_rc_file_is_read_in(self):
        assert "SIM_VERBOSE" in rc.__rc__

    def test_has_version_info(self):
        assert __version__


class TestRcFile:
    def test_rcfile_exists(self):
        assert os.path.exists(os.path.join(rc.__pkg_dir__, ".scopesimrc"))

    def test_rc_file_readable_by_scopesim_parser(self):
        rc_file =  os.path.join(rc.__pkg_dir__, ".scopesimrc")
        rc_dict = read_config(rc_file)
        assert isinstance(rc_dict, dict)
        assert len(rc_dict) > 0
