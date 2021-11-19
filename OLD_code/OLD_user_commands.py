import os
import logging
import copy

import yaml

from scopesim import rc
from scopesim.utils import find_file
from OLD_code import OLD_user_commands_utils as cutils

__all__ = ["UserCommands"]


class UserCommands:
    def __init__(self, filename=None, sim_data_dir=None, **kwargs):

        self.pkg_dir = rc.__pkg_dir__
        self.data_dir = rc.__data_dir__

        self._yaml_dicts = []

        # read in the default keywords
        # Don't use self.update because we need to add all the valid keywords
        self.cmds = copy.deepcopy(rc.__old_config__)
        self.cmds.update(rc.__rc__)
        self._convert_python_types()

        self.cmds["CONFIG_USER"] = filename
        self.cmds["CONFIG_DEFAULT"] = os.path.join(rc.__pkg_dir__,
                                                   "OLD.default.config")
        # read in the users wishes
        # Use self.update so that we reject all the invalid keywords
        if filename is not None:
            self.update(filename)

        # add any loose keywords
        self.update(kwargs)

        # option sim_data_dir overrides values in config files
        if sim_data_dir is not None:
            rc.__search_path__.insert(0, sim_data_dir)

    def update(self, new_input):
        if isinstance(new_input, UserCommands):
            tmp_cmds = new_input.cmds
        elif isinstance(new_input, dict):
            tmp_cmds = new_input
        elif isinstance(new_input, str):
            fname = find_file(new_input)
            if fname is not None:
                new_input = fname
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
        _yaml_dicts = []
        for yaml_name in ["SIM_GENERAL_YAML", "SIM_ATMOSPHERE_YAML",
                          "SIM_TELESCOPE_YAML", "SIM_RELAY_OPTICS_YAML",
                          "SIM_INSTRUMENT_YAML", "SIM_DETECTOR_YAML"]:
            yaml_obj = self.cmds[yaml_name]
            if isinstance(yaml_obj, dict):
                _yaml_dicts += [yaml_obj]
            elif isinstance(yaml_obj, str):
                with open(find_file(yaml_obj)) as f:
                    _yaml_dicts += [dic for dic in yaml.full_load_all(f)]

        self._yaml_dicts = _yaml_dicts
        return self._yaml_dicts

    def __getitem__(self, key):
        if cutils.is_item_subcategory(key, self.cmds):
            return cutils.get_subcategory(key, self.cmds)
        else:
            return self.cmds[key]

    def __setitem__(self, key, val):
        if key not in self.cmds:
            logging.warning("{} not in self.keys. Ignoring.".format(key))
            return None

        self.cmds[key] = cutils.str_to_python_type(val)
