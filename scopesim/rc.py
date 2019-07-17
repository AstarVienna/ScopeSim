import os
import yaml

from .commands.system_dict import SystemDict

__pkg_dir__ = os.path.dirname(__file__)
__data_dir__ = os.path.join(__pkg_dir__, "data")

with open(os.path.join(__pkg_dir__, "defaults.yaml")) as f:
    dicts = [dic for dic in yaml.load_all(f)]
__config__ = SystemDict(dicts)

__currsys__ = __config__
__search_path__ = ['./', __config__["!SIM.file.local_packages_path"],
                   __pkg_dir__, __data_dir__]   # For utils.find_file()
