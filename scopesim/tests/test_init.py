import sys
import os
import logging
import yaml

from scopesim import __version__, rc


class TestBasicLoading:

    def test_pkg_dir_in_search_path(self):
        assert rc.__pkg_dir__ in rc.__search_path__

    def test_defaults_config_file_is_read_in(self):
        assert rc.__config__["!SIM.reports.verbose"] is False

    def test_has_version_info(self):
        assert __version__


class TestDefaultsYamlFile:
    def test_default_yaml_file_exists(self):
        assert os.path.exists(os.path.join(rc.__pkg_dir__, "defaults.yaml"))

    def test_rc_file_readable_by_scopesim_parser(self):
        default_file = os.path.join(rc.__pkg_dir__, "defaults.yaml")
        with open(default_file, "r") as f:
            default_dict = [dic for dic in yaml.full_load_all(f)]

        assert isinstance(default_dict[0], dict)
        assert len(default_dict[0]) > 0
