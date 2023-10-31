# -*- coding: utf-8 -*-
"""Global configurations for ScopeSim."""

from pathlib import Path
import yaml

from copy import deepcopy

from .system_dict import SystemDict, UniqueList

__pkg_dir__ = Path(__file__).parent

with (__pkg_dir__ / "defaults.yaml").open(encoding="utf-8") as file:
    dicts = list(yaml.full_load_all(file))


try:
    with (Path.home() / ".scopesim_rc.yaml").open(encoding="utf-8") as file:
        dicts.extend(list(yaml.full_load_all(file)))
except FileNotFoundError:
    pass


__config__ = SystemDict(dicts)
__currsys__ = deepcopy(__config__)

# Order matters!
__search_path__ = UniqueList([
    Path(__config__["!SIM.file.local_packages_path"]).absolute(),
    Path(__pkg_dir__).absolute(),
    *[Path(pth).absolute() for pth in __config__["!SIM.file.search_path"]],
])
