from .. import effects as efs
from ..effects.effects_utils import make_effect, get_all_effects
from ..utils import clean_dict


class OpticalElement:
    def __init__(self, yaml_dict=None, **kwargs):
        self.meta = {"name": "<empty>"}
        self.meta.update(kwargs)
        self.properties = {}
        self.effects = []

        if isinstance(yaml_dict, dict):
            self.meta = {key: yaml_dict[key] for key in yaml_dict
                         if key not in ["properties", "effects"]}
            if "properties" in yaml_dict:
                self.properties = yaml_dict["properties"]
                self.properties = clean_dict(self.properties["kwargs"],
                                             self.meta)
            if "effects" in yaml_dict:
                self.effects_dicts = yaml_dict["effects"]
                self.make_effects(yaml_dict["effects"])

    def make_effects(self, effects_dicts):
        for effdic in effects_dicts:
            if "kwargs" in effdic:
                effdic["kwargs"] = clean_dict(effdic["kwargs"], self.meta)
            self.effects += [make_effect(effdic, **self.properties)]

    def add_effect(self, effect):
        if isinstance(effect, efs.Effect):
            self.effects += [effect]

    def get_all(self, effect_class):
        return get_all_effects(self.effects, effect_class)

    def get_z_order_effects(self, z_level):
        if isinstance(z_level, int):
            zmin = z_level
            zmax = zmin + 99
        elif isinstance(z_level, (tuple, list)):
            zmin, zmax = z_level[:2]
        else:
            zmin, zmax = 0, 600

        effects = []
        for eff in self.effects:
            if "z_order" in eff.meta:
                z = eff.meta["z_order"]
                if isinstance(z, (list, tuple)):
                    if any([zmin <= zi <= zmax for zi in z]):
                        effects += [eff]
                else:
                    if zmin <= z <= zmax:
                        effects += [eff]

        return effects

    @property
    def ter_list(self):
        ter_list = [effect for effect in self.effects
                    if isinstance(effect, (efs.SurfaceList, efs.TERCurve))]
        return ter_list

    @property
    def mask_list(self):
        mask_list = [effect for effect in self.effects
                     if isinstance(effect, efs.ApertureList)]
        return mask_list

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        return self.get_all(item)

    def __repr__(self):
        msg = '\nOpticalElement : "{}" contains {} Effects: \n' \
              ''.format(self.meta["name"], len(self.effects))
        eff_str = "\n".join(["[{}] {}".format(i, eff.__repr__())
                             for i, eff in enumerate(self.effects)])
        return msg + eff_str
