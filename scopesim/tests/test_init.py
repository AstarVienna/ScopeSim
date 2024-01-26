from pathlib import Path
import yaml

from scopesim import __version__, rc


DEAULT_YAML = Path(rc.__pkg_dir__, "defaults.yaml")


class TestBasicLoading:

    def test_pkg_dir_in_search_path(self):
        assert rc.__pkg_dir__ in rc.__search_path__

    def test_defaults_config_file_is_read_in(self):
        assert rc.__config__["!SIM.reports.verbose"] is False

    def test_has_version_info(self):
        assert __version__


class TestDefaultsYamlFile:
    def test_default_yaml_file_exists(self):
        assert DEAULT_YAML.exists()

    def test_rc_file_readable_by_scopesim_parser(self):
        with DEAULT_YAML.open("r", encoding="utf-8") as file:
            default_dict = list(yaml.full_load_all(file))

        assert isinstance(default_dict[0], dict)
        assert len(default_dict[0]) > 0
