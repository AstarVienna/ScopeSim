from ...source.source2 import Source
from ..fov import FieldOfView
from ..data_container import DataContainer


class Effect(DataContainer):

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = []

    def apply_to(self, obj):
        if obj is None:
            raise ValueError("fov was None. Oops.")

        if not isinstance(obj, (Source, FieldOfView)):
            raise ValueError("fov must be a FieldOfView (or Source) object: {}"
                             "".format(type(obj)))

        return obj

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"edges": None, "wavelengths": None}

    def update(self, **kwargs):
        pass

    def __repr__(self):
        if "name" not in self.meta:
            self.meta["name"] = "<unknown name>"
        return '{}: "{}"'.format(type(self).__name__, self.meta["name"])
