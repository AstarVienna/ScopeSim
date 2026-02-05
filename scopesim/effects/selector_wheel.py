"""
A new effect module that enables multi-arm/branch instruments by allowing application of different effects of same type
to different FoV objects. For e.g. for a simultaneous 2-arm instrument (red and blue), if aperture_id=1 corresponds to
FoV objects of blue arm and aperture_id=2 to FoV objects of the red arm (and they already have non-overlapping wavelengths),
most of the existing effects in the optical train drop the FoV objects that are outside their spatio-spectral volume.
There is a need of an effect that stores different effects of the same type that apply to different FoV objects based on
a selector (e.g. aperture_id).

This module implements a SelectorWheel effect that allows the user to define multiple effects of the same type
( e.g. different aperture masks for different arms) in the wheel dictionary where each effect corresponds to a
"selector_id" value. The user can set which id to use as the "selector", for e.g. aperture_id or id of the FoV object.
"""
from ..utils import (check_keys, get_logger)
from .effects import Effect
from ..optics.fov_volume_list import FovVolumeList
from ..optics.fov import FieldOfView
import importlib

logger = get_logger(__name__)


class SelectorWheel(Effect):
    """
    SelectorWheel class - that applies different effect objects (of the same type) based on a selector value.

    Examples:
    ---------
    ::
        name: aperture_selector_wheel
        class: SelectorWheel
        kwargs:
            selector_key: aperture_id
            wheel:
                - selector_value: 0
                  effect_class: ApertureList
                  effect_kwargs:
                      filename: "slits/slit_nir.txt"
                - selector_value: 1
                  effect_class: ApertureList
                  effect_kwargs:
                      filename: "slits/slit_vis.txt"

    In the above example, the SelectorWheel effect applies different aperture masks to FoV objects
    based on their aperture_id value. If a FoV object has aperture_id=0, it gets the NIR slit mask applied,
    while aperture_id=1 gets the VIS slit mask. This allows for multi-arm instruments to have different
    aperture masks applied in a single optical train. The effect_class specified in the wheel entries
    must be a valid Effect subclass available in scopesim.effects module.

    """

    z_order = (290, 690)
    required_keys = {"selector_key", "wheel"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_keys(kwargs, self.required_keys, action="error")

        self.meta.update(kwargs)

        self.wheel_effects = {}
        for wheel_entry in self.meta["wheel"]:
            selector_value = wheel_entry["selector_value"]
            effect_class_name = wheel_entry["effect_class"]
            effect_kwargs = wheel_entry.get("effect_kwargs", {})

            # Dynamically import the effect class
            effect_module = importlib.import_module("scopesim.effects")
            effect_class = getattr(effect_module, effect_class_name)
            # Instantiate the effect and store it in the wheel_effects dictionary
            self.wheel_effects[selector_value] = effect_class(**effect_kwargs)


    def apply_to(self, obj, **kwargs):
        """Based on the selector_key's value in obj.meta, apply the corresponding effect from the wheel."""
        if isinstance(obj, FieldOfView):
            if self.meta['selector_key'] not in obj.meta.keys():
                raise ValueError(f"Selector key {self.meta['selector_key']} not found in FieldOfView meta.")
            selector_value = obj.meta[self.meta["selector_key"]]

            effect_to_apply = self.get_effect(selector_value)
            if effect_to_apply is None:
                logger.warning(f"No effect found for selector value: {selector_value}, skipping effect application.")
                return obj

            return effect_to_apply.apply_to(obj, **kwargs)

        if isinstance(obj, FovVolumeList):
            unique_selector_values = set([vol["meta"].get(self.meta["selector_key"], None) for vol in obj.volumes])
            new_volumes = []

            for val in unique_selector_values:
                vols_with_val = [vol for vol in obj.volumes if vol["meta"].get(self.meta["selector_key"], None) == val]
                effect_to_apply = self.get_effect(val)
                print(f"Applying effect for selector value: {val} -> {effect_to_apply}, volumes: {len(vols_with_val)}")

                if val is None or effect_to_apply is None:
                    new_volumes.extend(vols_with_val)
                    continue

                #TODO: If effect_to_apply is a dichroic which reassigns selector_key values, we need to add a check here
                #TODO: i.e. dichroic.arm_action.keys() should not include items in unique_selector_values other than val

                newvollist = FovVolumeList()
                newvollist.volumes = vols_with_val
                newvollist = effect_to_apply.apply_to(newvollist, **kwargs)
                new_volumes.extend(newvollist.volumes)

            obj.volumes = new_volumes
            return obj


    def get_effect(self, selector_value):
        eff = None
        if (selector_value is None) or (selector_value not in self.wheel_effects.keys()):
            logger.warning(f"Either None or Missing value of {self.meta["selector_key"]} requested: {selector_value}, returning None.")
        else:
            eff = self.wheel_effects[selector_value]
        return eff


