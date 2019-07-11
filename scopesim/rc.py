import os
import yaml

from .commands.user_commands_utils import read_config
from .commands.system_dict import SystemDict

__pkg_dir__ = os.path.dirname(__file__)
__data_dir__ = os.path.join(__pkg_dir__, "data")

__system__ = SystemDict()

with open(os.path.join(__pkg_dir__, "defaults.yaml")) as f:
    for dic in yaml.load_all(f):
        __system__.update(dic)

# load in settings from config files
__rc__ = read_config(os.path.join(__pkg_dir__, ".scopesimrc"))
__config__ = read_config(os.path.join(__pkg_dir__, ".default.config"))

__search_path__ = ['./', __rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                   __pkg_dir__, __data_dir__]   # For utils.find_file()
