import os
import warnings
import copy

import numpy as np
import yaml

from .. import rc
from ..utils import find_file

__all__ = ["UserCommands"]


class UserCommands:
    """
    Contains all the setting that a user may wish to alter for an optical train

    Most of the important settings are kept in the ``.cmds`` system dictionary
    Setting can be accessed by using the alias names. Currently these are:

    - ATMO: atmospheric and observatory location settings
    - TEL: telescope related settings
    - RO: relay optics settings, i.e. between telescope and instrument
    - INST: instrument optics settings
    - DET: detector settings
    - OBS: observation settings, and
    - SIM: simulation settings

    All of the settings are contained in a special ``SystemDict`` dictionary
    that allows the user to access all the settings via a bang-string (!). E.g::

        cmds = UserCommands()
        cmds["!SIM.file.local_packages_path]

    .. note::
       To use this format for accessing hierarchically-stored values, the bang
       string must always begin with a "!"

    Alternatively the same value can be accessed via the normal dictionary
    format. E.g::

        cmds["SIM"]["file"]["local_packages_path"]


    Parameters
    ----------
    filename : str, optional
        The path to a YAML file containing a list of settings

    kwargs
    ------
    use_instrument : str, optional
        The name of the main instrument to use

    packages : list, optional
        list of package names needed for the optical system, so that ScopeSim
        can find the relevant files. E.g. ["Armazones", "ELT", "MICADO"]

    yamls : list, optional
        list of yaml filenames that are needed for the combined optical system
        E.g. ["MICADO_Standalone_RO.yaml", "MICADO_H4RG.yaml", "MICADO_.yaml"]

    properties : dict, optional
        Any extra "OBS" properties that should be added

    ignore_effects : list
        Not yet implemented

    add_effects : list
        Not yet implemented

    override_effect_values : dict
        Not yet implemented

    Attributes
    ----------
    cmds : SystemDict
        Built from the ``properties`` dictionary of a yaml dictionary. All
        values here are accessible globally by all ``Effects`` objects in an
        ``OpticalTrain`` once the ``UserCommands`` has been passed to the
        ``OpticalTrain``.

    yaml_dicts : list of dicts
        Where all the effects dictionaries are stored


    """

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
        """
        Updates the current parameters with a yaml dictionary

        Parameters
        ----------
        yaml_input : str, dict
            - str: path to a file containing a YAML dictionary
            - dict: a yaml style dict following the ``OpticalElement`` format

        """

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
        """
        Deals with all the extra keys that are allowed for an ``OBS`` yaml dict
        """

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
            if not os.path.exists(filename):
                raise ValueError("{} does not contain a default.yaml file. "
                                 "Please specify a package with a default.yaml"
                                 "".format(obs_dic["use_instrument"]))
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
    """
    Adds the paths of a list of locally saved packages to the search path list

    Parameters
    ----------
    local_path : str
        Where the pacakges are located. Generally given by the value in
        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    package_list : list
        A list of the package names to add

    """

    for pkg in package_list:
        pkg_dir = os.path.abspath(os.path.join(local_path, pkg))
        if not os.path.exists(pkg_dir):
            # todo: keep here, but add test for this by downloading test_package
            # raise ValueError("Package could not be found: {}".format(pkg_dir))
            warnings.warn("Package could not be found: {}".format(pkg_dir))
        rc.__search_path__ += [pkg_dir]


def load_yaml_dicts(filename):
    """
    Loads a the one or more dicts stored in a YAML file under ``filename``

    Parameters
    ----------
    filename : str
        Path to the YAML file

    Returns
    -------
    yaml_dicts : list
        A list of dicts

    """

    yaml_dicts = []
    with open(filename) as f:
        yaml_dicts += [dic for dic in yaml.load_all(f)]

    return yaml_dicts


def list_local_packages(action="display"):
    """
    Lists the packages on the local disk that ScopeSim can find

    Packages can only be found in the directory listed under::

        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    Packages are divided into "main" packages and "extension" packages.

    - Main packages contain a ``default.yaml`` file which tell ScopeSim which
      other packages are required to generate the full optical system
    - Extension packages contain only the data files needed to support the
      effects listed in the package YAML file

    .. note::
       Only "main" packages can be passed to a UserCommands object using the
       ``use_instrument=...`` parameter

    Parameters
    ----------
    action : str, optional
        ["display", "return"] What to do with the output.
        - "display": the list of packages are printed to the screen
        - "return": package names are returned in lists

    Returns
    -------
    main_pkgs, ext_pkgs : lists
        If action="return": Lists containing the names of locally saved packages

    """

    local_path = os.path.abspath(rc.__config__["!SIM.file.local_packages_path"])
    pkgs = [d for d in os.listdir(local_path) if
            os.path.isdir(os.path.join(local_path, d))]

    main_pkgs = [pkg for pkg in pkgs if
                 os.path.exists(os.path.join(local_path, pkg, "default.yaml"))]
    ext_pkgs = [pkg for pkg in pkgs if not
                os.path.exists(os.path.join(local_path, pkg, "default.yaml"))]

    if action == "display":
        msg = "\nLocal package directory:\n  {}\n" \
              "Full packages [can be used with 'use_instrument=...']\n  {}\n" \
              "Support packages\n  {}".format(local_path, main_pkgs, ext_pkgs)
        print(msg)
    else:
        return main_pkgs, ext_pkgs