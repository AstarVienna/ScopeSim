import os
import inspect
from .commands.user_commands_utils import read_config

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")
__search_path__ = ['./', __pkg_dir__, __data_dir__]   # For utils.find_file()

# load in settings from config files
__rc__ = read_config(os.path.join(__pkg_dir__, ".simscoperc"))
__config__ = read_config(os.path.join(__pkg_dir__, ".default.config"))
