import os
import yaml

from .system_dict import SystemDict

__pkg_dir__ = os.path.dirname(__file__)

with open(os.path.join(__pkg_dir__, "defaults.yaml")) as f:
    dicts = [dic for dic in yaml.full_load_all(f)]

user_rc_path = os.path.expanduser("~/.scopesim_rc.yaml")
if os.path.exists(user_rc_path):
    with open(user_rc_path) as f:
        dicts += [dic for dic in yaml.full_load_all(f)]

__config__ = SystemDict(dicts)
__currsys__ = __config__

__search_path__ = [__config__["!SIM.file.local_packages_path"],
                   __pkg_dir__] + __config__["!SIM.file.search_path"]

# if os.environ.get("READTHEDOCS") == "True" or "F:" in os.getcwd():
#     extra_paths = ["../", "../../", "../../../", "../../../../"]
#     __search_path__ = extra_pa    ths + __search_path__