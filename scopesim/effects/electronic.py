import numpy as np

from astropy.io import fits

from .. import rc
from . import Effect
from ..base_classes import DetectorBase
from ..utils import real_colname


class DarkCurrent(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [501]

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            if isinstance(self.meta["value"], dict):
                dtcr_id = obj.meta[real_colname("id", obj.meta)]
                dark = self.meta["value"][dtcr_id]
            elif isinstance(self.meta["value"], float):
                dark = self.meta["value"]
            else:
                raise ValueError("<DarkCurrent>.meta['value'] must be either"
                                 "dict or float: {}".format(self.meta["value"]))

            dit = obj.meta["OBS_DIT"]
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

            orig_type = type(obj.image_hdu.data[0, 0])
            if not isinstance(obj.image_hdu.data[0, 0], np.float64):
                obj.image_hdu.data = obj.image_hdu.data.astype(np.float64)

            data = obj.image_hdu.data
            # ..todo FIX THIS!!!!!!
            # data = np.random.poisson(data).astype(orig_type)

            new_imagehdu = fits.ImageHDU(data=data, header=obj.image_hdu.header)
            obj.image_hdu = new_imagehdu

        return obj


class LinearityCurve(Effect):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)


class ReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(ReadoutNoise, self).__init__(**kwargs)
