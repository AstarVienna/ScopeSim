from ..effects.data_container import DataContainer
from ..base_classes import SourceBase, FieldOfViewBase, ImagePlaneBase, \
    DetectorBase


class Effect(DataContainer):

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = []

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, (SourceBase, FieldOfViewBase,
                                ImagePlaneBase, DetectorBase)):
            raise ValueError("object must one of the following: "
                             "Source, FieldOfView, ImagePlane, Detector: "
                             "{}".format(type(obj)))

        return obj

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"edges": None, "wavelengths": None}

    def update(self, **kwargs):
        pass

    def __repr__(self):
        if "name" not in self.meta:
            self.meta["name"] = "<unknown name>"
        return '{}: "{}"'.format(type(self).__name__, self.meta["name"])
