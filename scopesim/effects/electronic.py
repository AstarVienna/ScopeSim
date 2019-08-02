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
        self.meta["random_seed"] = rc.__currsys__["!SIM.random.seed"]
        self.meta.update(kwargs)

        self.required_keys = ["noise_std", "n_channels", "ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det):
        if isinstance(det, DetectorBase):
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"])
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            from_currsys(self.meta)
            ron_keys = ["noise_std", "n_channels", "channel_fraction",
                        "line_fraction", "pedestal_fraction", "read_fraction"]
            ron_kwargs = {key: self.meta[key] for key in ron_keys}
            ron_kwargs["image_shape"] = det.image_hdu.data.shape

            for n in range(self.meta["ndit"]):
                det.image_hdu.data += make_ron_frame(**ron_kwargs)

        return det


def make_ron_frame(image_shape, noise_std, n_channels, channel_fraction,
                   line_fraction, pedestal_fraction, read_fraction):
    shape = image_shape
    w_chan = shape[0] // n_channels

    pixel_std = noise_std * (pedestal_fraction + read_fraction)**0.5
    line_std = noise_std * line_fraction**0.5
    if shape < (1024, 1024):
        pixel = np.random.normal(loc=0, scale=pixel_std, size=shape)
        line = np.random.normal(loc=0, scale=line_std, size=shape[1])
    else:
        pixel = pseudo_random_field(scale=pixel_std, size=shape)
        line = pixel[0]

    channel_std = noise_std * channel_fraction**0.5
    channel = np.repeat(np.random.normal(loc=0, scale=channel_std,
                                         size=n_channels),
                        w_chan, axis=0)

    ron_frame = (pixel + line).T + channel

    return ron_frame


def pseudo_random_field(scale=1, size=(1024, 1024)):
    n = 256
    image = np.zeros(size)
    batch = np.random.normal(loc=0, scale=scale, size=(2*n, 2*n))
    for y in range(0, size[1], n):
        for x in range(0, size[0], n):
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
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"])
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

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
            ndit = from_currsys(self.meta["ndit"])

            obj.image_hdu.data += dark * dit * ndit

        return obj


class LinearityCurve(Effect):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [540]

        self.required_keys = ["ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det):
        if isinstance(det, DetectorBase):
            ndit = from_currsys(self.meta["ndit"])
            incident = self.table["incident"] * ndit
            measured = self.table["measured"] * ndit

            image = det.image_hdu.data
            shape = image.shape
            flat_image = image.flatten()
            new_flat_image = np.interp(flat_image, incident, measured)
            new_image = new_flat_image.reshape(shape)

            det.image_hdu.data = new_image

        return det




