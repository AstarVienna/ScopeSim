import numpy as np

from astropy.io import fits

from .. import rc
from . import Effect
from ..base_classes import DetectorBase, ImagePlaneBase
from ..utils import real_colname, from_currsys, check_keys, interp2


class SummedExposure(Effect):
    def __init__(self, **kwargs):
        super(SummedExposure, self).__init__(**kwargs)
        self.meta["z_order"] = [860]

        required_keys = ["dit", "ndit"]
        check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj):
        if isinstance(obj, DetectorBase):
            dit = from_currsys(self.meta["dit"])
            ndit = from_currsys(self.meta["ndit"])

            obj._hdu.data *= dit * ndit

        return obj


class PoorMansHxRGReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(PoorMansHxRGReadoutNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [811]
        self.meta["pedestal_fraction"] = 0.3
        self.meta["read_fraction"] = 0.4
        self.meta["line_fraction"] = 0.25
        self.meta["channel_fraction"] = 0.05
        self.meta["random_seed"] = "!SIM.random.seed"
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
            ron_kwargs["image_shape"] = det._hdu.data.shape

            for _ in range(self.meta["ndit"]):
                det._hdu.data += make_ron_frame(**ron_kwargs)

        return det


class BasicReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(BasicReadoutNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [811]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

        self.required_keys = ["noise_std", "ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det):
        if isinstance(det, DetectorBase):
            self.meta = from_currsys(self.meta)
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            noise_std = self.meta["noise_std"] * np.sqrt(float(self.meta["ndit"]))
            det._hdu.data += np.random.normal(loc=0, scale=noise_std,
                                              size=det._hdu.data.shape)

        return det


class ShotNoise(Effect):
    def __init__(self, **kwargs):
        super(ShotNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [820]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

    def apply_to(self, det):
        if isinstance(det, DetectorBase):
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"])
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            # ! poisson(x) === normal(mu=x, sigma=x**0.5)
            # Windows has a porblem with generating poisson values above 2**30
            # Above ~100 counts the poisson and normal distribution are
            # basically the same. For large arrays the normal distribution
            # takes only 60% as long as the poisson distribution
            data = det._hdu.data
            below = data < 2**20
            above = np.invert(below)
            data[below] = np.random.poisson(data[below]).astype(float)
            data[above] = np.random.normal(data[above], np.sqrt(data[above]))
            data = np.floor(data)
            new_imagehdu = fits.ImageHDU(data=data, header=det._hdu.header)
            det._hdu = new_imagehdu

        return det


class DarkCurrent(Effect):
    """
    required: dit, ndit, value
    """
    def __init__(self, **kwargs):
        super(DarkCurrent, self).__init__(**kwargs)
        self.meta["z_order"] = [830]

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

            obj._hdu.data += dark * dit * ndit

        return obj


class LinearityCurve(Effect):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)
        self.meta["z_order"] = [840]

        self.required_keys = ["ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det):
        if isinstance(det, DetectorBase):
            ndit = from_currsys(self.meta["ndit"])
            incident = self.table["incident"] * ndit
            measured = self.table["measured"] * ndit

            image = det._hdu.data
            shape = image.shape
            flat_image = image.flatten()
            new_flat_image = np.interp(flat_image, incident, measured)
            new_image = new_flat_image.reshape(shape)

            det._hdu.data = new_image

        return det


class ReferencePixelBorder(Effect):
    def __init__(self, **kwargs):
        super(ReferencePixelBorder, self).__init__(**kwargs)
        self.meta["z_order"] = [780]
        val = 0
        if "all" in kwargs:
            val = int(kwargs["all"])
        widths = {key: val for key in ["top", "bottom", "left", "right"]}
        self.meta.update(widths)
        self.meta.update(kwargs)

    def apply_to(self, implane):
        if isinstance(implane, ImagePlaneBase):
            if self.meta["top"] > 0:
                implane.hdu.data[:, -self.meta["top"]:] = 0
            if self.meta["bottom"] > 0:
                implane.hdu.data[:, :self.meta["bottom"]] = 0
            if self.meta["right"] > 0:
                implane.hdu.data[-self.meta["right"]:, :] = 0
            if self.meta["left"] > 0:
                implane.hdu.data[:self.meta["left"], :] = 0

        return implane


################################################################################


def make_ron_frame(image_shape, noise_std, n_channels, channel_fraction,
                   line_fraction, pedestal_fraction, read_fraction):
    shape = image_shape
    w_chan = max(1, shape[0] // n_channels)

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
                                         size=n_channels), w_chan + 1, axis=0)

    ron_frame = (pixel + line).T + channel[:shape[0]]

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
