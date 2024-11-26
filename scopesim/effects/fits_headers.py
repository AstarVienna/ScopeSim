# -*- coding: utf-8 -*-
"""TBA."""

from copy import deepcopy
import datetime
from typing import ClassVar

import yaml
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from astar_utils.nested_mapping import recursive_update

from . import Effect
from ..source.source_fields import HDUSourceField, TableSourceField
from ..utils import from_currsys, find_file


class ExtraFitsKeywords(Effect):
    """
    Extra FITS header keywords to be added to the pipeline FITS files.

    These keywords are ONLY for keywords that should be MANUALLY ADDED to the
    headers after a simulation is read-out by the detector.

    Simulation parameters (Effect kwargs values, etc) will be added
    automatically by ScopeSim in a different function, but following this
    format.

    The dictionaries should be split into different HIERARCH lists, e.g.:

    - HIERARCH ESO
      For ESO specific keywords
    - HIERARCH SIM
      For ScopeSim specific keywords, like simulation parameters
    - HIERARCH MIC
      For MICADO specific keywords, (unsure what these would be yet)

    More HIERARCH style keywords can also be added as needed for other
    use-cases.

    Parameters
    ----------
    filename : str, optional
        Name of a .yaml nested dictionary file. See below for examples

    yaml_string : str, optional
        A triple-" string containing the contents of a yaml file

    header_dict : nested dicts, optional
        A series of nested python dictionaries following the format of the
        examples below. This keyword allows these dicts to be definied directly
        in the Effect yaml file, rather than in a seperate header keywords
        file.

    Examples
    --------
    Specifying the extra FITS keywords directly in the .yaml file where the
    Effect objects are described.
    ::

        name: extra_fits_header_entries
        class: ExtraFitsKeywords
        kwargs:
          header_dict:
            - ext_type: PrimaryHDU
              keywords:
                HIERARCH:
                  ESO:
                    ATM:
                      TEMPERAT: -5

    The contents of ``header_dict`` can also be abstracted away into a seperate
    file, e.g. ``extra_FITS_keys.yaml``. The file format is described below in
    detail below.
    ::

        name: extra_fits_header_entries
        class: ExtraFitsKeywords
        kwargs:
          filename: extra_FITS_keys.yaml

    The Effect can be added directly in an iPython session.
    ::

        >>> hdr_dic = {"ext_type": "PrimaryHDU",
                       "keywords":
                           {"HIERARCH":
                               {"SIM":
                                   {"hello": world}
                               }
                           }
                       }
        >>> extra_keys = ExtraFitsKeywords(header_dict=hdr_dic)
        >>> optical_train.optics_manager.add_effect(extra_keys)


    Yaml file format
    ----------------
    This document is a yaml document.
    Hence all new keywords should be specified in the form of yaml nested
    dictionaries.
    As each ``astropy.HDUList`` contains one or more extensions, the inital
    level is reserved for a list of keyword groups.
    For example::

        - ext_type: PrimaryHDU
          keywords:
            HIERARCH:
              ESO:
                ATM:
                  TEMPERAT: -5

        - ext_number: [1, 2]
          keywords:
            HIERARCH:
              ESO:
                DET:
                  DIT: [5, '[s] exposure length']   # example of adding a comment
            EXTNAME: "DET§.DATA"                   # example of extension specific qualifier

    The keywords can be added to one or more extensions, based on one of the
    following ``ext_`` qualifiers: ``ext_name``, ``ext_number``, ``ext_type``

    Each of these ``ext_`` qualifiers can be a ``str`` or a ``list``.
    For a list, ScopeSim will add the keywords to all extensions matching the
    specified type/name/number

    The number of the extension can be used in a value by using the "§" string.
    That is, keyword values with "§" with have the extension number inserted
    where the "§" is.

    The above example (``EXTNAME: "DET§.DATA"``) will result in the following
    keyword added only to extensions 1 and 2:

    - PrimaryHDU (ext 0)::

          header['HIERARCH ESO ATM TEMPERAT'] = -5

    - Extension 1 (regardless of type)::

          header['HIERARCH ESO DET DIT'] = (5, '[s] exposure length')
          header['EXTNAME'] = "DET1.DATA"

    - Extension 2 (regardless of type)::

          header['HIERARCH ESO DET DIT'] = (5, '[s] exposure length')
          header['EXTNAME'] = "DET2.DATA"

    Resolved and un-resolved keywords
    ---------------------------------
    ScopeSim uses bang-strings to resolve global parameters.
    E.g: ``from_currsys('!ATMO.temperature')`` will resolve to a float
    These bang-strings will be resolved automatically in the ``keywords``
    dictionary section.

    If the keywords bang-string should instead remain unresolved and the string
    added verbatim to the header, we use the ``unresolved_keywords`` dictionary
    section.

    Additionally, new functionality will be added to ScopeSim to resolve the
    kwargs/meta parameters of Effect objects.
    The format for this will be to use a new type: the hash-string.
    This will have this format::

        #<optical_element_name>.<effect_name>.<kwarg_name>

    For example, the temperature of the MICADO detector array can be accessed
    by::

        '#MICADO_DET.full_detector_array.temperature'

    In the context of the yaml file this would look like::

        - ext_type: PrimaryHDU
          keywords:
            HIERARCH:
              ESO:
                DET
                  TEMPERAT: '#MICADO_DET.full_detector_array.temperature'

    Obviously some though needs to be put into how exactly we list the
    simulation parameters in a coherent manner.
    But this is 'Zukunftsmusik'.
    For now we really just want an interface that can add the ESO header
    keywords, which can also be expanded in the future for our own purposes.

    Below is an example of some extra keywords for MICADO headers::

        - ext_type: PrimaryHDU
          keywords:
            HIERARCH:
              ESO:
                ATM:
                  TEMPERAT: '!ATMO.temperature'   # will be resolved via from_currsys
                  PWV: '!ATMO.pwv'
                  SEEING: 1.2
                DAR:
                  VALUE: '#<effect_name>.<meta_name>'   # will be resolved via effects
                DPR:
                  TYPE: 'some_type'
              SIM:
                random_simulation_keyword: some_value
              MIC:
                micado_specific: ['keyword', 'keyword comment']

          unresolved_keywords:
            HIERARCH:
              ESO:
                ATM:
                  TEMPERAT: '!ATMO.temperature'   # will be left as a string

        - ext_type: ImageHDU
          keywords:
            HIERARCH:
              SIM:
                hello: world
                hallo: welt
                grias_di: woed
                zdrasviute: mir
                salud: el mundo

    """

    z_order: ClassVar[tuple[int, ...]] = (1050,)

    def __init__(self, cmds=None, **kwargs):
        # don't pass kwargs, as DataContainer can't handle yaml files
        super().__init__(cmds=cmds)
        params = {"name": "extra_fits_keywords",
                  "description": "Extra FITS headers",
                  "header_dict": None,
                  "filename": None,
                  "yaml_string": None,
                  }
        self.meta.update(params)
        self.meta.update(kwargs)

        tmp_dicts = []
        if self.meta["filename"] is not None:
            yaml_file = find_file(self.meta["filename"])
            with open(yaml_file, encoding="utf-8") as file:
                # possible multiple yaml docs in a file
                # --> returns list even for a single doc
                tmp_dicts.extend(dic for dic in yaml.full_load_all(file))

        if self.meta["yaml_string"] is not None:
            yml = self.meta["yaml_string"]
            tmp_dicts.extend(dic for dic in yaml.full_load_all(yml))

        if self.meta["header_dict"] is not None:
            if not isinstance(self.meta["header_dict"], list):
                tmp_dicts.append(self.meta["header_dict"])
            else:
                tmp_dicts.extend(self.meta["header_dict"])

        self.dict_list = []
        for dic in tmp_dicts:
            # format says yaml file contains list of dicts
            if isinstance(dic, list):
                self.dict_list.extend(dic)
            # catch case where user forgets the list
            elif isinstance(dic, dict):
                self.dict_list.append(dic)

    def apply_to(self, hdul, **kwargs):
        """
        Add extra fits keywords from a yaml file including !,#-stings.

        Parameters
        ----------
        optical_train : scopesim.OpticalTrain, optional
            Used to resolve #-strings

        """
        opt_train = kwargs.get("optical_train")
        if isinstance(hdul, fits.HDUList):
            for dic in self.dict_list:
                resolved = flatten_dict(dic.get("keywords", {}), resolve=True,
                                        optics_manager=opt_train)
                unresolved = flatten_dict(dic.get("unresolved_keywords", {}))
                exts = get_relevant_extensions(dic, hdul)
                for i in exts:
                    # On windows machines Â appears in the string when using §
                    resolved_with_counters = {
                        k: v.replace("Â",
                                     "").replace("§",
                                                 str(i)).replace("++", str(i))
                        if isinstance(v, str) else v for k, v in resolved.items()
                    }
                    hdul[i].header.update(resolved_with_counters)
                    hdul[i].header.update(unresolved)

        return hdul


def get_relevant_extensions(dic, hdul):
    exts = []
    if dic.get("ext_name") is not None:
        exts.extend(i for i, hdu in enumerate(hdul)
                    if hdu.header["EXTNAME"] == dic["ext_name"])
    elif dic.get("ext_number") is not None:
        ext_n = np.array(dic["ext_number"])
        exts.extend(ext_n[ext_n < len(hdul)])
    elif dic.get("ext_type") is not None:
        if isinstance(dic["ext_type"], list):
            ext_type_list = dic["ext_type"]
        else:
            ext_type_list = [dic["ext_type"]]
        cls = tuple(getattr(fits, cls_str) for cls_str in ext_type_list)
        exts.extend(i for i, hdu in enumerate(hdul) if isinstance(hdu, cls))

    return exts


def flatten_dict(dic, base_key="", flat_dict=None, resolve=False,
                 optics_manager=None, cmds=None):
    """
    Flattens nested yaml dictionaries into a single level dictionary.

    Parameters
    ----------
    dic : dict
    base_key : str
    flat_dict : dict, optional
        Top-level dictionary for recursive calls
    resolve : bool
        If True, resolves !-str via from_currsys and #-str via optics_manager
    optics_manager : scopesim.OpticsManager
        Required for resolving #-strings
    cmds : UserCommands
        To use for resolving !-strings

    Returns
    -------
    flat_dict : dict

    """
    if cmds is None and optics_manager is not None:
        cmds = optics_manager.cmds

    if flat_dict is None:
        flat_dict = {}
    for key, val in dic.items():
        flat_key = f"{base_key}{key} "
        if isinstance(val, dict):
            flatten_dict(val, flat_key, flat_dict, resolve, optics_manager, cmds)
        else:
            flat_key = flat_key[:-1]

            # catch any value+comments lists
            comment = ""
            if isinstance(val, list) and len(val) == 2 and isinstance(val[1],
                                                                      str):
                value, comment = val
            else:
                value = deepcopy(val)

            # resolve any bang or hash strings
            if resolve and isinstance(value, str):
                if value.startswith("!"):
                    value = from_currsys(value, cmds)
                elif value.startswith("#"):
                    if optics_manager is None:
                        raise ValueError("An OpticsManager object must be "
                                         "passed in order to resolve "
                                         "#-strings")
                    value = optics_manager[value]

            if isinstance(value, u.Quantity):
                comment = f"[{str(value.unit)}] {comment}"
                value = value.value

            # Convert e.g.  Unit(mag) to just "mag". Not sure how this will
            # work when deserializing though.
            if isinstance(value, u.Unit):
                value = str(value)

            if isinstance(value, (list, np.ndarray)):
                value = f"{value.__class__.__name__}:{str(list(value))}"
                max_len = 80 - len(flat_key)
                if len(value) > max_len:
                    value = f"{value[:max_len-4]} ..."

            if isinstance(value, (datetime.time, datetime.date,
                                  datetime.datetime)):
                value = value.isoformat()

            # Add the flattened KEYWORD = (value, comment) to the header dict
            if comment:
                flat_dict[flat_key] = (value, str(comment))
            else:
                flat_dict[flat_key] = value

    return flat_dict


class EffectsMetaKeywords(ExtraFitsKeywords):
    """
    Adds meta dictionary info from all Effects to the FITS headers.

    Parameters
    ----------
    ext_number : int, list of ints, optional
        Default 0. The numbers of the extensions to which the header keywords
        should be added

    add_excluded_effects : bool, optional
        Default False. Add meta dict for effects with
        ``<effect>.include=False``

    keyword_prefix : str, optional
        Default "HIERARCH SIM". Custom FITS header keyword prefix. Effect meta
        dict entries will appear in the header as:
        ``<keyword_prefix> EFFn <key> : <value>``

    Examples
    --------
    Yaml file entry:
    ::

        name: effect_dumper
        class: EffectsMetaKeywords
        description: adds all effects meta dict entries to the FITS header
        kwargs:
          ext_number: [0, 1]
          add_excluded_effects: False
          keyword_prefix: HIERARCH SIM

    """

    z_order: ClassVar[tuple[int, ...]] = (1040,)

    def __init__(self, cmds=None, **kwargs):
        super(ExtraFitsKeywords, self).__init__(cmds=cmds, **kwargs)
        params = {"name": "effects_fits_keywords",
                  "description": "Effect Meta FITS headers",
                  "ext_number": [0],
                  "add_excluded_effects": False,
                  "keyword_prefix": "HIERARCH SIM"}
        self.meta.update(params)
        self.meta.update(kwargs)

    def apply_to(self, hdul, **kwargs):
        """See parent docstring."""
        opt_train = kwargs.get("optical_train")
        if isinstance(hdul, fits.HDUList) and opt_train is not None:
            # todo: use a different way of getting all the effect names
            # opt.effects returns the __repr__, not the original name
            for i, eff_name in enumerate(opt_train.effects["name"]):
                # Check for spaces
                if " " in eff_name:
                    # E.g. 'filter_wheel_1 : [open]'
                    assert "wheel" in eff_name, \
                        f"Unknown effect name with space: {eff_name}"
                    eff_name = eff_name.split()[0]

                # get a resolved meta dict from the effect
                eff_meta = deepcopy(opt_train[eff_name].meta)

                if self.meta["add_excluded_effects"] and not eff_meta["include"]:
                    continue

                keys = list(eff_meta.keys())
                for key in keys:
                    value = eff_meta[key]
                    if key in ["history", "notes", "changes", "cmds"]:
                        eff_meta.pop(key)
                    if isinstance(value, Table):
                        eff_meta[key] = f"Table object of length: {len(value)}"

                # add effect under the EFFn keyword
                prefix = self.meta["keyword_prefix"]
                class_name = opt_train[eff_name].__class__.__name__
                self.dict_list = [
                    {"ext_number": self.meta["ext_number"],
                     "keywords": {
                         f"{prefix} EFF{i} class": [
                             class_name,
                             "ScopeSim class name"
                             ],
                         f"{prefix} EFF{i}": eff_meta
                         }
                     }
                    ]
                hdul = super().apply_to(hdul=hdul, optical_train=opt_train)

        return hdul


class SourceDescriptionFitsKeywords(ExtraFitsKeywords):
    """
    Adds parameters from all Source fields to the FITS headers.

    Parameters
    ----------
    ext_number : int, list of ints, optional
        Default 0. The numbers of the extensions to which the header keywords
        should be added

    keyword_prefix : str, optional
        Default "HIERARCH SIM". Custom FITS header keyword prefix. Effect meta
        dict entries will appear in the header as:
        ``<keyword_prefix> SRCn <key> : <value>``

    Examples
    --------
    Yaml file entry:
    ::

        name: source_descriptor
        class: SourceDescriptionFitsKeywords
        description: adds info from all Source fields to the FITS header
        kwargs:
          ext_number: [0]
          keyword_prefix: HIERARCH SIM

    """

    z_order: ClassVar[tuple[int, ...]] = (1030,)

    def __init__(self,  cmds=None, **kwargs):
        super(ExtraFitsKeywords, self).__init__(cmds=cmds, **kwargs)
        params = {"name": "source_fits_keywords",
                  "description": "Source description FITS headers",
                  "ext_number": [0],
                  "keyword_prefix": "HIERARCH SIM"}
        self.meta.update(params)
        self.meta.update(kwargs)

    def apply_to(self, hdul, **kwargs):
        """See parent docstring."""
        opt_train = kwargs.get("optical_train")
        if not isinstance(hdul, fits.HDUList) or opt_train is None:
            return hdul

        if (src := opt_train._last_source) is not None:
            prefix = self.meta["keyword_prefix"]
            for i, field in enumerate(src.fields):
                src_class = field.field.__class__.__name__
                src_dic = deepcopy(src.meta)
                if isinstance(field, HDUSourceField):
                    hdr = field.header
                    for key in hdr:
                        src_dic = {key: [hdr[key], hdr.comments[key]]}

                elif isinstance(field, TableSourceField):
                    src_dic.update(field.meta)
                    src_dic["length"] = len(field)
                    for j, name in enumerate(field.field.colnames):
                        src_dic[f"col{j}_name"] = name
                        src_dic[f"col{j}_unit"] = str(field.field[name].unit)

                self.dict_list = [{"ext_number": self.meta["ext_number"],
                                   "keywords": {
                                       f"{prefix} SRC{i} class": src_class,
                                       f"{prefix} SRC{i}": src_dic}
                                   }]
                hdul = super().apply_to(hdul=hdul, optical_train=opt_train)

        # catch the function call
        for hdu in hdul:
            for key in hdu.header:
                if "function_call" in key:
                    hdu.header[f"FN{key.split()[1]}"] = hdu.header.pop(key)

        return hdul


class SimulationConfigFitsKeywords(ExtraFitsKeywords):
    """
    Adds parameters from all config dictionaries to the FITS headers.

    Parameters
    ----------
    ext_number : int, list of ints, optional
        Default 0. The numbers of the extensions to which the header keywords
        should be added

    resolve : bool
        Default True. If True, all !-strings and #-strings are resolved via
        ``from_currsys`` before being add to the header. If False, the
        unaltered !-strings or #-strings are added to the header.

    keyword_prefix : str, optional
        Default "HIERARCH SIM". Custom FITS header keyword prefix. Effect meta
        dict entries will appear in the header as:
        ``<keyword_prefix> SRCn <key> : <value>``

    Examples
    --------
    Yaml file entry:
    ::

        name: source_descriptor
        class: SimulationConfigFitsKeywords
        description: adds info from all config dicts to the FITS header
        kwargs:
          ext_number: [0]
          resolve: False
          keyword_prefix: HIERARCH SIM

    """

    z_order: ClassVar[tuple[int, ...]] = (1020,)

    def __init__(self, cmds=None, **kwargs):
        super(ExtraFitsKeywords, self).__init__(cmds=cmds, **kwargs)
        params = {"name": "simulation_fits_keywords",
                  "description": "Simulation Config FITS headers",
                  "ext_number": [0],
                  "resolve": True,
                  "keyword_prefix": "HIERARCH SIM"}
        self.meta.update(params)
        self.meta.update(kwargs)

    def apply_to(self, hdul, **kwargs):
        """See parent docstring."""
        opt_train = kwargs.get("optical_train")
        if isinstance(hdul, fits.HDUList) and opt_train is not None:
            # HACK: This workaround was added after the ChainMap change, to
            #       have a simply but save way of getting the final dict out.
            # TODO: Improve this at some point.....
            cmds = deepcopy(opt_train.cmds.maps[-1].dic)
            for m in opt_train.cmds.maps[-2::-1]:
                cmds = recursive_update(cmds, m.dic)

            sim_prefix = self.meta["keyword_prefix"]
            resolve_prefix = "unresolved_" if not self.meta["resolve"] else ""
            # needed for the super().apply_to method
            self.dict_list = [{"ext_number": self.meta["ext_number"],
                               f"{resolve_prefix}keywords": {
                                   f"{sim_prefix} CONFIG": cmds}
                               }]
            hdul = super().apply_to(hdul=hdul, optical_train=opt_train)

        return hdul
