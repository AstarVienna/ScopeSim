import numpy as np

from .effects import Effect
from ..base_classes import SourceBase, ImagePlaneBase, FieldOfViewBase
from ..utils import required_keys, check_keys

class RotateSource(Effect):
    def __init__(self, **kwargs):
        super(Rotation2D, self).__init__(**kwargs)
        params = {"z_order": [530],
                  "report_plot_include": False,
                  "report_table_include": False, }
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ["angle"]
        check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, SourceBase):
            obj.rotate(angle=self.meta["angle"])

        return obj

