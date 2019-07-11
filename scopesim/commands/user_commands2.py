import os
import warnings
import copy

import yaml

from .. import rc
from ..utils import find_file
from . import user_commands_utils as cutils

__all__ = ["UserCommands"]


class UserCommands:
    def __init__(self, filename=None, sim_data_dir=None, **kwargs):

        self.filename = filename
        self.meta = {}
        self.meta.update(kwargs)
        self._yaml_dicts = []
        rc.SystemDict["SIM.file.search_path"] += ["sim_data_dir"]

        # read in the default keywords
        # Don't use self.update because we need to add all the valid keywords

    @property
    def yaml_dicts(self):
        _yaml_dicts = []
        for yaml_name in ["SIM_GENERAL_YAML", "SIM_ATMOSPHERE_YAML",
                          "SIM_TELESCOPE_YAML", "SIM_RELAY_OPTICS_YAML",
                          "SIM_INSTRUMENT_YAML", "SIM_DETECTOR_YAML"]:
            yaml_obj = self.cmds[yaml_name]
            if isinstance(yaml_obj, dict):
                _yaml_dicts += [yaml_obj]
            elif isinstance(yaml_obj, str):
                with open(find_file(yaml_obj)) as f:
                    _yaml_dicts += [dic for dic in yaml.load_all(f)]

        self._yaml_dicts = _yaml_dicts
        return self._yaml_dicts

