import os
import yaml

from .commands.user_commands_utils import read_config
from .commands.system_dict import SystemDict

__pkg_dir__ = os.path.dirname(__file__)
__data_dir__ = os.path.join(__pkg_dir__, "data")

with open(os.path.join(__pkg_dir__, "defaults.yaml")) as f:
    dicts = [dic for dic in yaml.load_all(f)]
__config__ = SystemDict(dicts)
__currsys__ = __config__

# load in settings from config files
__rc__ = read_config(os.path.join(__pkg_dir__, ".scopesimrc"))
__old_config__ = read_config(os.path.join(__pkg_dir__, ".default.config"))

__search_path__ = ['./', __rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                   __pkg_dir__, __data_dir__]   # For utils.find_file()
