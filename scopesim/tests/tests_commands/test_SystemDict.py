import pytest
import copy
import yaml

from scopesim.system_dict import SystemDict, recursive_update

_basic_yaml = """
alias : OBS
properties :
    temperature : 100    
"""


@pytest.fixture(scope="class")
def basic_yaml():
    return yaml.full_load(_basic_yaml)


@pytest.mark.usefixtures("basic_yaml")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(SystemDict(), SystemDict)

    def test_initalises_with_normal_dict(self):
        sys_dict = SystemDict({"a": 1})
        assert isinstance(sys_dict, SystemDict)

    def test_initalises_with_yaml_dict(self, basic_yaml):
        sys_dict = SystemDict(basic_yaml)
        assert isinstance(sys_dict, SystemDict)
        assert "OBS" in sys_dict.dic


@pytest.mark.usefixtures("basic_yaml")
class TestActsLikeDict:
    def test_can_add_and_retrieve_normal_dict_entries(self):
        sys_dict = SystemDict()
        sys_dict["name"] = "ELT"
        assert sys_dict["name"] == "ELT"

    def test_can_and_and_retrieve_special_dict_entries(self, basic_yaml):
        sys_dict = SystemDict(basic_yaml)
        sys_dict["!OBS.lam.max.unit"] = "um"
        print(sys_dict)
        assert sys_dict["!OBS.lam.max.unit"] == "um"
        assert sys_dict["!OBS.temperature"] == 100

    def test_uses___contains___keyword_for_normal_dicts(self):
        sys_dict = SystemDict({"name": "ELT"})
        assert "name" in sys_dict

    def test_uses___contains___keyword_for_special_dicts(self, basic_yaml):
        sys_dict = SystemDict(basic_yaml)
        assert "!OBS.temperature" in sys_dict
        assert "!OBS.temperature.unit" not in sys_dict
        assert "!OBS.humidity" not in sys_dict


@pytest.mark.usefixtures("basic_yaml")
class TestRecursiveUpdate:
    def test_updates_normal_recursive_dicts(self):
        sys_dict = SystemDict()
        sys_dict["name"] = {"value": "ELT"}
        sys_dict.update({"name": {"type": "str"}})
        assert sys_dict["name"]["value"] == "ELT"
        assert sys_dict["name"]["type"] == "str"

    def test_updates_yaml_alias_recursive_dicts(self, basic_yaml):
        sys_dict = SystemDict(copy.deepcopy(basic_yaml))
        basic_yaml["properties"] = {"temperature": 42, "humidity": 0.75}
        sys_dict.update(basic_yaml)
        assert sys_dict["!OBS.temperature"] == 42
        assert sys_dict["OBS"]["humidity"] == 0.75


class TestFunctionRecursiveUpdate:
    def test_recursive_update_combines_dicts(self):
        e = {"a": {"b": {"c": 1}}}
        f = {"a": {"b": {"d": 2}}}
        recursive_update(e, f)
        assert e["a"]["b"]["c"] == 1
        assert e["a"]["b"]["d"] == 2

    def test_recursive_update_overwrites_scalars(self):
        e = {"a": {"b": {"c": 1}}}
        f = {"a": {"b": {"c": 2}}}
        recursive_update(e, f)
        assert e["a"]["b"]["c"] == 2

    def test_recursive_update_overwrites_dict_with_scalar(self):
        e = {"a": {"b": {"c": 1}}}
        f = {"a": {"b": 4}}
        recursive_update(e, f)
        assert e["a"]["b"] == 4

    def test_recursive_update_overwrites_scalar_with_dict(self):
        e = {"a": {"b": 5}}
        f = {"a": {"b": {"c": 1}}}
        recursive_update(e, f)
        assert e["a"]["b"] == {"c": 1}

    def test_recursive_update_overwrites_string_with_string(self):
        e = {"a": {"b": {"c": "hello"}}}
        f = {"a": {"b": {"c": "world"}}}
        recursive_update(e, f)
        assert e["a"]["b"]["c"] == "world"
