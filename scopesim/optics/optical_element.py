import logging
from inspect import isclass
from typing import TextIO
from io import StringIO

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
                              if key not in {"properties", "effects"}})
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

                    self.effects.append(make_effect(eff_dic, **self.properties))

    def add_effect(self, effect):
        if isinstance(effect, efs.Effect):
            self.effects.append(effect)
        else:
            logging.warning("%s is not an Effect object and was not added", effect)

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
                    if any(zmin <= zi <= zmax for zi in z):
                        effects.append(eff)
                else:
                    if zmin <= z <= zmax:
                        effects.append(eff)

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
        """
        Returns Effects of Effect meta properties

        Parameters
        ----------
        item : str, int, Effect-Class
            Either the name, list index, or Class of an effect.
            It is possible to get meta entries from a named effect by using:
            ``#<effect_name>.<meta-key>``

        Returns
        -------
        obj : str, float, Effect, list
            - str, float : if #-string is used to get Effect.meta entries
            - Effect : if a unique Effect name is given
            - list : if an Effect-Class is given

        Examples
        --------
        ::
            from scopesim.effect import TERCurve
            opt_el = opt_elem.OpticalElement(detector_yaml_dict)

            effect_list = opt_el[TERCurve]
            effect = opt_el[0]
            effect = opt_el["detector_qe_curve"]
            meta_value = opt_el["#detector_qe_curve.filename"]

        """
        obj = None
        if isclass(item):
            obj = self.get_all(item)
        elif isinstance(item, int):
            obj = self.effects[item]
        elif isinstance(item, str):
            if item.startswith("#") and "." in item:
                eff, meta = item.replace("#", "").split(".")
                obj = self[eff][f"#{meta}"]
            else:
                obj = [eff for eff in self.effects
                       if eff.meta["name"] == item]

        if isinstance(obj, list) and len(obj) == 1:
            obj = obj[0]
        # if obj is None or len(obj) == 0:
        #     logging.warning(f'No result for key: "{item}". '
        #                     f'Did you mean "#{item}"?')

        return obj

    def write_string(self, stream: TextIO, list_effects: bool = True) -> None:
        """Write formatted string representation to I/O stream"""
        stream.write(f"{self!s} contains {len(self.effects)} Effects\n")
        if list_effects:
            for i_eff, eff in enumerate(self.effects):
                stream.write(f"[{i_eff}] {eff!r}\n")

    def pretty_str(self) -> str:
        """Return formatted string representation as str"""
        with StringIO() as str_stream:
            self.write_string(str_stream)
            output = str_stream.getvalue()
        return output

    @property
    def display_name(self):
        return self.meta.get("name", self.meta.get("filename", "<empty>"))

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return f"{self.__class__.__name__}: \"{self.display_name}\""

    @property
    def properties_str(self):
        prop_str = ""
        max_key_len = max(len(key) for key in self.properties.keys())
        padlen = max_key_len + 4
        for key in self.properties:
            if key not in {"comments", "changes", "description", "history",
                           "report"}:
                prop_str += f"{key:>{padlen}} : {self.properties[key]}\n"

        return prop_str

    def report(self, filename=None, output="rst", rst_title_chars="^#*+",
               **kwargs):

        rst_str = f"""
{str(self)}
{rst_title_chars[0] * len(str(self))}

**Element**: {self.meta.get("object", "<unknown optical element>")}

**Alias**: {self.meta.get("alias", "<unknown alias>")}
        
**Description**: {self.meta.get("description", "<no description>")}

Global properties
{rst_title_chars[1] * 17}
::

{self.properties_str}
"""

        if len(self.list_effects()) > 0:
            rst_str += f"""
Effects
{rst_title_chars[1] * 7}

Summary of Effects included in this optical element:

.. table::
    :name: {"tbl:" + self.meta.get("name", "<unknown OpticalElement>")}
   
{table_to_rst(self.list_effects(), indent=4)}
 
"""

        reports = [eff.report(rst_title_chars=rst_title_chars[-2:], **kwargs)
                   for eff in self.effects]
        rst_str += "\n\n" + "\n\n".join(reports)

        write_report(rst_str, filename, output)

        return rst_str
