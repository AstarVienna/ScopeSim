from copy import deepcopy, copy

from astropy.table import Table

from .. import effects as efs


def combine_surface_effects_OLD(surface_effects):
    surflist_list = [eff for eff in surface_effects
                     if isinstance(eff, efs.SurfaceList)]
    surf_list = [eff for eff in surface_effects
                 if isinstance(eff, efs.TERCurve)]

    if len(surflist_list) == 0:
        surf = empty_surface_list()
        surf.meta["name"] = "Radiometry Table"
        surflist_list += [surf]

    new_surflist = deepcopy(surflist_list[0])
    for surflist in surflist_list[1:]:
        new_surflist.add_surface_list(surflist)

    for surf in surf_list:
        position = surf.meta["position"] if "position" in surf.meta else -1
        new_surflist.add_surface(surf, surf.meta["name"], position=position)

    new_surflist.table = new_surflist.radiometry_table.table

    return new_surflist


def combine_surface_effects(surface_effects):
    surflist_list = [eff for eff in surface_effects
                     if isinstance(eff, efs.SurfaceList)]
    surf_list = [eff for eff in surface_effects
                 if isinstance(eff, (efs.TERCurve, efs.FilterWheel))
                 and not isinstance(eff, efs.SurfaceList)]

    if len(surflist_list) == 0:
        surflist_list = [empty_surface_list(name="combined_surface_list")]

    new_surflist = copy(surflist_list[0])
    for surflist in surflist_list[1:]:
        new_surflist.add_surface_list(surflist)

    # ..todo:: should read position from the list positions in surface_effects
    for surf in surf_list:
        position = surf.meta.get("position", -1)
        new_surflist.add_surface(surf, surf.meta["name"], position=position)

    return new_surflist


def get_all_effects(effects, effect_class):
    if isinstance(effect_class, (list, tuple)):
        my_effects = []
        for eff_cls in effect_class:
            my_effects += get_all_effects(effects, eff_cls)
    else:
        my_effects = [eff for eff in effects
                      if isinstance(eff, effect_class) and eff.include]

    return my_effects


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
    return any([isinstance(eff, efs.SpectralTraceList) for eff in effects])


def empty_surface_list(**kwargs):
    tbl = Table(names=["name", "outer", "inner", "angle",
                       "temperature", "action", "filename"],
                data=[["test"], [0.], [0.], [0.], [0.], ["none"], ["none"]],
                meta={"outer_unit": "m", "inner_unit": "m",
                      "angle_unit": "deg", "temperature_unit": "deg_C"})
    return efs.SurfaceList(table=tbl[:0], **kwargs)
