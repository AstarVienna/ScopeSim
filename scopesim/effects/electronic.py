import numpy as np

from astropy.io import fits

from .. import rc
from . import Effect
from ..base_classes import DetectorBase
from ..utils import real_colname, from_currsys, check_keys


class BasicReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(BasicReadoutNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [510]
        self.meta["pedestal_fraction"] = 0.3
        self.meta["read_fraction"] = 0.4
        self.meta["line_fraction"] = 0.25
        self.meta["channel_fraction"] = 0.05
        self.meta.update(kwargs)

        self.required_keys = ["noise_std", "n_channels"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det):
        if isinstance(det, DetectorBase):
            from_currsys(self.meta)
            ron_keys = ["ndit", "noise_std", "n_channels", "channel_fraction",
                        "line_fraction", "pedestal_fraction", "read_fraction"]
            ron_kwargs = {self.meta[key] for key in ron_keys}
            ron_kwargs["image_shape"] = det.image_hdu.data.shape

            det.image_hdu.data += make_ron_frame(**self.meta)

        return det


def make_ron_frame(image_shape, ndit, noise_std, n_channels, channel_fraction,
                   line_fraction, pedestal_fraction, read_fraction):
    shape = list(image_shape) + [ndit]
    w_chan = shape[0] // n_channels
    pixel = np.random.random(shape) * (pedestal_fraction + read_fraction)
    line = np.random.random(shape[1:]) * line_fraction
    channel = np.repeat(np.random.random([n_channels, ndit]),
                        w_chan, axis=0) * channel_fraction

    total = (np.rot90(pixel + line, axes=(0, 1)) + channel) * noise_std
    total_flat = np.sum(total, axis=2)

    return total_flat


def pseudo_random_field(shape=(1024, 1024)):
    n = 256
    image = np.zeros(shape)
    batch = np.random.random((2*n, 2*n))
    for y in range(0, shape[1], n):
        for x in range(0, shape[0], n):
            i, j = np.random.randint(n, size=2)
            image[x:x+n, y:y+n] = batch[i:i+n, j:j+n]

    return image


class ShotNoise(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [520]
        self.meta["random_seed"] = rc.__currsys__["!SIM.random.seed"]
        self.meta.update(kwargs)

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            seed = self.meta["random_seed"]
            if isinstance(seed, str) and seed[0] == "!":
                seed = from_currsys(seed)
                self.meta["random_seed"] = seed

            if seed is not None:
                np.random.seed(seed)

            orig_type = type(det.image_hdu.data[0, 0])
            if not isinstance(det.image_hdu.data[0, 0], np.float64):
                det.image_hdu.data = det.image_hdu.data.astype(np.float64)

            data = det.image_hdu.data
            # ..todo FIX THIS!!!
            # (KL) it seems to have fixed itself... lets wait and see
            data = np.random.poisson(data).astype(orig_type)

            new_imagehdu = fits.ImageHDU(data=data, header=det.image_hdu.header)
            det.image_hdu = new_imagehdu

        return det


class DarkCurrent(Effect):
    """
    required: dit, ndit, value

    """
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [530]

        required_keys = ["value", "dit", "ndit"]
        check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj):
        if isinstance(obj, DetectorBase):
            if isinstance(self.meta["value"], dict):
                dtcr_id = obj.meta[real_colname("id", obj.meta)]
                dark = self.meta["value"][dtcr_id]
            elif isinstance(self.meta["value"], float):
                dark = self.meta["value"]
            else:
                raise ValueError("<DarkCurrent>.meta['value'] must be either"
                                 "dict or float: {}".format(self.meta["value"]))

            dit = from_currsys(self.meta["dit"])
            obj.image_hdu.data += dark * dit

        return obj


class LinearityCurve(Effect):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [540]
