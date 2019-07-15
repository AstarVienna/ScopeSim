import os
import warnings
import copy

import numpy as np
import yaml

from .system_dict import SystemDict
from .. import rc
from ..utils import find_file

__all__ = ["UserCommands"]


class UserCommands:
    """

    kwargs
    ------
    name:



    """

    def __init__(self, filename=None, **kwargs):

        self.cmds = copy.deepcopy(rc.SystemDict)
        self.yaml_dicts = []

        self.filename = filename
        if filename is not None:
            self.yaml_dicts += load_yaml_dicts(find_file(filename))
        else:
            kwargs["alias"] = "OBS"
            self.yaml_dicts += [kwargs]

        self.update()

    def update(self, **kwargs):
        # combine dicts with same alias
        # update self.cmds
        self.yaml_dicts = combine_similar_yamls(self.yaml_dicts)



        for yaml_dic in self.yaml_dicts:
            self.cmds.update(yaml_dic)



def combine_similar_yamls(yaml_list):
    aliases = np.unique([yaml["alias"] for yaml in yaml_list])
    for alias in alias




    return yaml_list



def load_yaml_dicts(filename):
    yaml_dicts = []
    with open(filename) as f:
        yaml_dicts += [dic for dic in yaml.load_all(f)]

    return yaml_dicts

