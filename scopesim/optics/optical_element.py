import logging
from inspect import isclass

from astropy.table import Table

from .. import effects as efs
from ..effects.effects_utils import make_effect, get_all_effects
from ..utils import write_report
from ..reports.rst_utils import table_to_rst
from .. import rc


class OpticalElement:
    """
    Contains all information to describe a section of an optical system

    There are 5 major section: ``location``, ``telescope``, ``relay optics``
    ``instrument``, ``detector``.

    An OpticalElement describes how a certain section of the optical train
    changes the incoming photon distribution by specifying a list of
    ``Effects`` along with a set of local properties e.g. temperature, etc
    which are common to more than one Effect


    Parameters
    ----------
    yaml_dict : dict
        Description of optical section properties, effects, and meta-data

    kwargs : dict
        Optical Element specific information which has no connection to the
        effects that are passed. Any global values, e.g. airmass
        (i.e. bang strings) are passed on to the individual effect, where the
        relevant values are then pulled from rc.__currsys__


    Attributes
    ----------
    meta : dict
        Contains meta data from the yaml, and is updated with the OBS_DICT
        key-value pairs

    properties : dict
        Contains any properties that is "global" to the optical element, e.g.
        instrument temperature, atmospheric pressure. Any OBS_DICT keywords are
        cleaned with from the ``meta`` dict during initialisation

    effects : list of dicts
        Contains the a list of dict descriptions of the effects that the optical
        element generates. Any OBS_DICT keywords are cleaned with from the
        ``meta`` dict during initialisation


    Examples
    --------

    """
    def __init__(self, yaml_dict=None, **kwargs):
        self.meta = {"name": "<empty>"}
        self.meta.update(kwargs)
        self.properties = {}
        self.effects = []

        if isinstance(yaml_dict, dict):
            self.meta.update({key: yaml_dict[key] for key in yaml_dict
                              if key not in ["properties", "effects"]})
            if "properties" in yaml_dict:
                self.properties = yaml_dict["properties"]
            if "name" in yaml_dict:
                self.properties["element_name"] = yaml_dict["name"]
            if "effects" in yaml_dict and len(yaml_dict["effects"]) > 0:
                for eff_dic in yaml_dict["effects"]:
                    if "name" in eff_dic and hasattr(rc.__currsys__,
                                                     "ignore_effects"):
                        if eff_dic["name"] in rc.__currsys__.ignore_effects:
                            eff_dic["include"] = False

                    self.effects += [make_effect(eff_dic, **self.properties)]

    def add_effect(self, effect):
        if isinstance(effect, efs.Effect):
            self.effects += [effect]
        else:
            logging.warning("{} is not an Effect object and was not added"
                          "".format(effect))

    def get_all(self, effect_class):
        return get_all_effects(self.effects, effect_class)

    def get_z_order_effects(self, z_level):
        if isinstance(z_level, int):
            zmin = z_level
            zmax = zmin + 99
        elif isinstance(z_level, (tuple, list)):
            zmin, zmax = z_level[:2]
        else:
            zmin, zmax = 0, 999

        effects = []
        for eff in self.effects:
            if eff.include and "z_order" in eff.meta:
                z = eff.meta["z_order"]
                if isinstance(z, (list, tuple)):
                    if any([zmin <= zi <= zmax for zi in z]):
                        effects += [eff]
                else:
                    if zmin <= z <= zmax:
                        effects += [eff]

        return effects

    @property
    def surfaces_list(self):
        _ter_list = [effect for effect in self.effects
                     if isinstance(effect, (efs.SurfaceList, efs.FilterWheel,
                                            efs.TERCurve))]
        return _ter_list

    @property
    def masks_list(self):
        _mask_list = [effect for effect in self.effects if
                      isinstance(effect, (efs.ApertureList, efs.ApertureMask))]
        return _mask_list

    def list_effects(self):
        elements = [self.meta["name"]] * len(self.effects)
        names = [eff.display_name for eff in self.effects]
        classes = [eff.__class__.__name__ for eff in self.effects]
        included = [eff.meta["include"] for eff in self.effects]
        z_orders = [eff.meta["z_order"] for eff in self.effects]

        colnames = ["element", "name", "class", "included", "z_orders"]
        data = [elements, names, classes, included, z_orders]
        tbl = Table(names=colnames, data=data, copy=False)

        return tbl

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        obj = None
        if isclass(item):
            obj = self.get_all(item)
        elif isinstance(item, int):
            obj = self.effects[item]
        elif isinstance(item, str):
            obj = [eff for eff in self.effects
                   if eff.meta["name"] == item]

        return obj

    def __repr__(self):
        msg = '\nOpticalElement : "{}" contains {} Effects: \n' \
              ''.format(self.meta["name"], len(self.effects))
        eff_str = "\n".join(["[{}] {}".format(i, eff.__repr__())
                             for i, eff in enumerate(self.effects)])
        return msg + eff_str

    def __str__(self):
        name = self.meta.get("name", self.meta.get("filename", "<empty>"))
        return '{}: "{}"'.format(type(self).__name__, name)

    @property
    def properties_str(self):
        prop_str = ""
        max_key_len = max([len(key) for key in self.properties.keys()])
        for key in self.properties:
            if key not in ["comments", "changes", "description", "history",
                           "report"]:
                prop_str += "    {} : {}\n".format(key.rjust(max_key_len),
                                                   self.properties[key])

        return prop_str

    def report(self, filename=None, output="rst", rst_title_chars="^#*+",
               **kwargs):

        rst_str = """
{}
{}

**Element**: {}

**Alias**: {}
        
**Description**: {}

Global properties
{}
::

{}
""".format(str(self),
           rst_title_chars[0] * len(str(self)),
           self.meta.get("object", "<unknown optical element>"),
           self.meta.get("alias", "<unknown alias>"),
           self.meta.get("description", "<no description>"),
           rst_title_chars[1] * 17,
           self.properties_str)

        if len(self.list_effects()) > 0:
            rst_str += """        
Effects
{}

Summary of Effects included in this optical element:

.. table::
    :name: {}
   
{}
 
""".format(rst_title_chars[1] * 7,
           "tbl:" + self.meta.get("name", "<unknown OpticalElement>"),
           table_to_rst(self.list_effects(), indent=4))

        reports = [eff.report(rst_title_chars=rst_title_chars[-2:], **kwargs)
                   for eff in self.effects]
        rst_str += "\n\n" + "\n\n".join(reports)

        write_report(rst_str, filename, output)

        return rst_str
