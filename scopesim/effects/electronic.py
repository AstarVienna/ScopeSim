from . import Effect
from ..detector import Detector


class DarkCurrent(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [501]

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, Detector):
            dark = self.meta["value"]
            dit = self.meta["OBS_DIT"]
            obj.image_hdu.data += dark * dit

        return obj
