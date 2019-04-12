import os
import glob
import warnings
import copy

import yaml
import numpy as np
import astropy.io.ascii as ioascii

from .. import rc
from ..utils import find_file
from .. import server as svr
from . import user_commands_utils as cutils

__all__ = ["UserCommands"]


class UserCommands:
    def __init__(self, filename=None, sim_data_dir=None,
                 instrument=None, mode=None, filter_name=None):

        self.pkg_dir = rc.__pkg_dir__
        self.data_dir = rc.__data_dir__

        # read in the default keywords
        # Don't use self.update because we need to add all the valid keywords
        self.cmds = copy.deepcopy(rc.__config__)
        self.cmds.update(rc.__rc__)
        self._convert_python_types()

        self.cmds["CONFIG_USER"] = filename
        self.cmds["CONFIG_DEFAULT"] = os.path.join(rc.__pkg_dir__,
                                                   ".default.config")
        # read in the users wishes
        # Use self.update so that we reject all the invalid keywords
        if filename is not None:
            self.update(filename)

        # option sim_data_dir overrides values in config files
        if sim_data_dir is not None:
            self.cmds['SIM_DATA_DIR'] = sim_data_dir

        if instrument is not None:
            self.set_instrument(instrument)

        if mode is not None:
            self.set_mode(mode)

        if filter_name is not None:
            self.select_filter(filter_name)

    def update(self, new_input):
        if isinstance(new_input, UserCommands):
            tmp_cmds = new_input.cmds
        elif isinstance(new_input, dict):
            tmp_cmds = new_input
        elif isinstance(new_input, str):
            tmp_cmds = cutils.read_config(new_input)
        else:
            raise ValueError("Cannot update with type: " + type(new_input))

        # Use self.update so that we reject all the invalid keywords
        for key in tmp_cmds:
            self[key] = tmp_cmds[key]

    def _convert_python_types(self):
        """Convert strings "none", "true", "false" values into python types"""
        self.cmds = cutils.convert_dict_strings_to_python_types(self.cmds)

    @property
    def yaml_dicts(self):
        atmo_yaml = find_file(self.cmds["SIM_ATMOSPHERE_YAML"])
        scope_yaml = find_file(self.cmds["SIM_TELESCOPE_YAML"])
        inst_yaml = find_file(self.cmds["SIM_INSTRUMENT_YAML"])
        detector_yaml = find_file(self.cmds["SIM_DETECTOR_YAML"])

        yaml_dicts = []
        for yaml_file in [atmo_yaml, scope_yaml, inst_yaml, detector_yaml]:
            if yaml_file is not None:
                with open(find_file(yaml_file)) as f:
                    yaml_dicts += [dic for dic in yaml.load_all(f)]

        return yaml_dicts

    def __getitem__(self, key):
        if cutils.is_item_subcategory(key, self.cmds):
            return cutils.get_subcategory(key, self.cmds)
        else:
            return self.cmds[key]

    def __setitem__(self, key, val):
        if key not in self.cmds:
            warnings.warn("{} not in self.keys. Ignoring.".format(key))
            return None

        self.cmds[key] = cutils.str_to_python_type(val)
