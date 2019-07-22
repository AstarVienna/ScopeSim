from copy import deepcopy

from .. import effects as efs
from ..optics.radiometry_utils import empty_surface_list


def combine_surface_effects(surface_effects):
    surflist_list = [eff for eff in surface_effects
                     if isinstance(eff, efs.SurfaceList)]
    surf_list = [eff for eff in surface_effects
                 if isinstance(eff, efs.TERCurve)]

    if len(surflist_list) == 0:
        tbl = empty_surface_list()
        tbl.meta["name"] = "Radiometry Table"
        surflist_list += [tbl]

    new_surflist = deepcopy(surflist_list[0])
    for surflist in surflist_list[1:]:
        new_surflist.add_surface_list(surflist)

    for surf in surf_list:
        new_surflist.add_surface(surf, surf.meta["name"])

    return new_surflist


def get_all_effects(effects, effect_class):
    return [eff for eff in effects if isinstance(eff, effect_class)]


def make_effect(effect_dict, **properties):
    effect_meta_dict = {key : effect_dict[key] for key in effect_dict
                        if key not in ["class", "kwargs"]}
    effect_class_name = effect_dict["class"]
    effect_cls = getattr(efs, effect_class_name)
    # ..todo: add looking for custom effect class names from 3rd party packages

    effect_kwargs = {}
    effect_kwargs.update(effect_meta_dict)         # effect name and description
    effect_kwargs.update(properties)              # optical_element properties
    if "kwargs" in effect_dict:
        effect_kwargs.update(effect_dict["kwargs"])  # individual effect kwargs

    effect = effect_cls(**effect_kwargs)
    # effect.meta.update(effect_meta_dict)  # is this needed? Seems redundant

    return effect


def is_spectroscope(effects):
    has_trace_lists = sum([isinstance(eff, efs.SpectralTraceList)
                           for eff in effects])
    has_apertures = sum([isinstance(eff, (efs.ApertureList, efs.ApertureMask))
                         for eff in effects])

    return bool(has_apertures and has_trace_lists)
