"""TBA."""

import inspect
from copy import deepcopy, copy
from collections.abc import Iterable

from astropy.table import Table

from .. import effects as efs
from ..utils import get_logger


logger = get_logger(__name__)


# TODO: is this ever used anywhere??
def combine_surface_effects(surface_effects):
    surflist_list = [eff for eff in surface_effects
                     if isinstance(eff, efs.SurfaceList)]
    surf_list = [eff for eff in surface_effects
                 if isinstance(eff, (efs.TERCurve, efs.FilterWheelBase))
                 and not isinstance(eff, efs.SurfaceList)]

    if not surflist_list:
        surflist_list = [empty_surface_list(name="combined_surface_list")]

    new_surflist = copy(surflist_list[0])
    new_surflist.data_container = copy(surflist_list[0].data_container)

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


def make_effect(effect_dict, cmds=None, **properties):
    effect_meta_dict = {key: effect_dict[key] for key in effect_dict
                        if key not in ["class", "kwargs"]}
    effect_class_name = effect_dict["class"]
    effect_cls = getattr(efs, effect_class_name)
    # ..todo: add looking for custom effect class names from 3rd party packages

    effect_kwargs = {}
    effect_kwargs.update(effect_meta_dict)        # effect name and description
    effect_kwargs.update(properties)              # optical_element properties
    if "kwargs" in effect_dict:
        effect_kwargs.update(effect_dict["kwargs"])  # individual effect kwargs

    effect = effect_cls(cmds=cmds, **effect_kwargs)

    return effect


def is_spectroscope(effects):
    spec_classes = (efs.SpectralTraceList,
                    efs.SpectralTraceListWheel)
    return any(isinstance(eff, spec_classes) for eff in effects)


def empty_surface_list(**kwargs):
    tbl = Table(names=["name", "outer", "inner", "angle",
                       "temperature", "action", "filename"],
                data=[["test"], [0.], [0.], [0.], [0.], ["none"], ["none"]],
                meta={"outer_unit": "m", "inner_unit": "m",
                      "angle_unit": "deg", "temperature_unit": "deg_C"})
    return efs.SurfaceList(table=tbl[:0], **kwargs)


def scopesim_effect_classes(base_effect=efs.Effect):
    members = inspect.getmembers(efs)
    efs_dict = {".".join([cls.__module__, cls.__name__]).replace("scopesim.effects.", ""): cls
                for name, cls in members
                if hasattr(cls, "__mro__") and base_effect in cls.__mro__}
    sorted_effects = {key: efs_dict[key] for key in sorted(efs_dict)}

    return sorted_effects


def z_order_in_range(z_eff, z_range: range) -> bool:
    """
    Return True if any of the z_orders in `z_eff` is in the given range.

    The `z_range` parameter can be constructed as ``range(z_min, z_max)``.

    Parameters
    ----------
    z_eff : int or list of ints
        z_order(s) of the effect.
    z_range : range
        range object of allowed z_order values.

    Returns
    -------
    bool
        True if at least one z_order is in range, False otherwise.

    """
    if not isinstance(z_eff, Iterable):
        logger.warning("z_order %d should be a single-item iterable", z_eff)
        z_eff = [z_eff]

    return any(zi in z_range for zi in z_eff)
