import os
import inspect
from .commands.user_commands_utils import read_config

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# load in settings from config files
__rc__ = read_config(os.path.join(__pkg_dir__, ".scopesimrc"))
__config__ = read_config(os.path.join(__pkg_dir__, ".default.config"))

__search_path__ = ['./', __rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                   __pkg_dir__, __data_dir__]   # For utils.find_file()


class SystemDict(object):
    def __init__(self, new_dict=None):
        self.dic = {}
        if new_dict is not None:
            self.update(new_dict)

    def __getitem__(self, item):
        if isinstance(item, str) and item[0] == "!":
            item_chunks = item[1:].split(".")
            entry = self.dic
            for item in item_chunks:
                entry = entry[item]
            return entry
        else:
            return self.dic[item]

    def __setitem__(self, key, value):
        if isinstance(key, str) and key[0] == "!":
            key_chunks = key[1:].split(".")
            entry = self.dic
            for key in key_chunks[:-1]:
                if key not in entry:
                    entry[key] = {}
                entry = entry[key]
            entry[key_chunks[-1]] = value
        else:
            self.dic[key] = value

    def update(self, new_dict):
        if isinstance(new_dict, dict) and "alias" in new_dict:
            self.dic[new_dict["alias"]] = new_dict["properties"]
        else:
            self.dic.update(new_dict)


def test_system_dict():

    import yaml
    new_dict = yaml.load("""
alias : OBS
properties :
    temperature : 100    
    """)

    system_dict = SystemDict(new_dict)
    assert system_dict["!OBS.temperature"] == 100

    system_dict["!OBS.lam.max.unit"] = "um"
    assert system_dict["!OBS.lam.max.unit"] == "um"

    system_dict["name"] = "ELT"
    assert system_dict["name"] == "ELT"

    print(system_dict.dic)
