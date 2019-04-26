import numpy as np

from .. import rc
from . import Effect
from ..base_classes import DetectorBase


class DarkCurrent(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [501]

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            if isinstance(self.meta["value"], dict):
                dtcr_id = obj.meta["id"]
                dark = self.meta["value"][dtcr_id]
            elif isinstance(self.meta["value"], float):
                dark = self.meta["value"]
            else:
                raise ValueError("<DarkCurrent>.meta['value'] must be either"
                                 "dict or float: {}".format(self.meta["value"]))

            dit = self.meta["OBS_DIT"]
            obj.image_hdu.data += dark * dit

        return obj


class ShotNoise(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [502]

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):

            if "random_seed" in kwargs:
                if kwargs["random_seed"] == "SIM_RANDOM_SEED":
                    seed = rc.__rc__["SIM_RANDOM_SEED"]
                elif isinstance(kwargs["random_seed"], float):
                    seed = kwargs["random_seed"]
                else:
                    raise ValueError("random_seed must be a float or the RC "
                                     "keyword 'SIM_RANDOM_SEED'")
                np.random.seed(seed)

            if not isinstance(obj.image_hdu.data, np.float64):
                obj.image_hdu.data = np.random.poisson(obj.image_hdu.data)
            else:
                orig_type = type(obj.image_hdu.data)
                obj.image_hdu.data = obj.image_hdu.data.astype(np.float64)
                obj.image_hdu.data = np.random.poisson(obj.image_hdu.data)
                obj.image_hdu.data = obj.image_hdu.data.astype(orig_type)

        return obj


class LinearityCurve(Effect):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)


class ReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(ReadoutNoise, self).__init__(**kwargs)
