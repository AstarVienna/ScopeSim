from pathlib import Path
import yaml

from .system_dict import SystemDict

__pkg_dir__ = Path(__file__).parent

with open(__pkg_dir__/"defaults.yaml") as f:
    dicts = list(yaml.full_load_all(f))

user_rc_path = Path("~/.scopesim_rc.yaml").expanduser()
if user_rc_path.exists():
    with open(user_rc_path) as f:
        dicts.extend(list(yaml.full_load_all(f)))

__config__ = SystemDict(dicts)
__currsys__ = __config__

__search_path__ = [__config__["!SIM.file.local_packages_path"],
                   __pkg_dir__] + __config__["!SIM.file.search_path"]

# if os.environ.get("READTHEDOCS") == "True" or "F:" in os.getcwd():
#     extra_paths = ["../", "../../", "../../../", "../../../../"]
#     __search_path__ = extra_paths + __search_path__
