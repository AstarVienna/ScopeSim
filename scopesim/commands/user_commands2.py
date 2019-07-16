import os
import warnings
import copy

import numpy as np
import yaml

from .. import rc
from ..utils import find_file

__all__ = ["UserCommands"]


class UserCommands:

    def __init__(self, filename=None, **kwargs):

        self.cmds = copy.deepcopy(rc.__config__)
        self.yaml_dicts = []
        self.filename = filename

        self.ignore_effects = {}

        _yaml_dicts = []
        if filename is not None:
            _yaml_dicts += load_yaml_dicts(find_file(filename))
        else:
            kwargs["alias"] = "OBS"
            _yaml_dicts += [kwargs]

        for _yaml_dict in _yaml_dicts:
            self.update(_yaml_dict)

    def update(self, yaml_input):

        if isinstance(yaml_input, str):
            yaml_dicts = load_yaml_dicts(find_file(yaml_input))
        elif isinstance(yaml_input, dict):
            yaml_dicts = [yaml_input]
        else:
            raise ValueError("yaml_dicts must be a filename or a dictionary: {}"
                             "".format(yaml_input))

        self.yaml_dicts += yaml_dicts
        for yaml_dic in yaml_dicts:
            self.cmds.update(yaml_dic)
            if "alias" in yaml_dic and yaml_dic["alias"] == "OBS":
                self._import_obs_dict(yaml_dic)

    def _import_obs_dict(self, obs_dic):

        local_path = self.cmds["!SIM.file.local_packages_path"]
        if "packages" in obs_dic:
            add_packages_to_rc_search(local_path, obs_dic["packages"])

        if "yamls" in obs_dic:
            yaml_dicts = []
            for filename in obs_dic["yamls"]:
                yaml_dicts += load_yaml_dicts(find_file(filename))
            for yaml_dict in yaml_dicts:
                self.update(yaml_dict)

        if "use_instrument" in obs_dic:
            filename = os.path.join(os.path.abspath(local_path),
                                    obs_dic["use_instrument"],
                                    "default.yaml")
            for yaml_dict in load_yaml_dicts(filename):
                self.update(yaml_dict)

        if "ignore_effects" in obs_dic:
            self.ignore_effects = obs_dic["ignore_effects"]

        if "add_effects" in obs_dic:
            # ..todo: implement this
            pass
        if "override_effect_values" in obs_dic:
            # ..todo: implement this
            pass

    def __setitem__(self, key, value):
        self.cmds.__setitem__(key, value)

    def __getitem__(self, item):
        return self.cmds.__getitem__(item)

    def __contains__(self, item):
        return self.cmds.__contains__(item)

    def __repr__(self):
        return self.cmds.__repr__()


def add_packages_to_rc_search(local_path, package_list):

    for pkg in package_list:
        pkg_dir = os.path.abspath(os.path.join(local_path, pkg))
        if not os.path.exists(pkg_dir):
            # todo: keep here, but add test for this by downloading test_package
            # raise ValueError("Package could not be found: {}".format(pkg_dir))
            warnings.warn("Package could not be found: {}".format(pkg_dir))
        rc.__search_path__ += [pkg_dir]


def load_yaml_dicts(filename):

    yaml_dicts = []
    with open(filename) as f:
        yaml_dicts += [dic for dic in yaml.load_all(f)]

    return yaml_dicts

