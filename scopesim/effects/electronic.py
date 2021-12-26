"""
Electronic detector effects - related to detector readout

Classes:
- AutoExposure - determine DIT and NDIT automatically
- SummedExposure - simulates a summed stack of ``ndit`` exposures
- PoorMansHxRGReadoutNoise - simple readout noise for HAWAII detectors
- BasicReadoutNoise - readout noise
- ShotNoise - realisation of Poissonian photon noise
- DarkCurrent - add dark current
- LinearityCurve - apply detector (non-)linearity and saturation
- ReferencePixelBorder
- BinnedImage

Functions:
- make_ron_frame
- pseudo_random_field
"""
import numpy as np

from astropy.io import fits

from .. import rc
from . import Effect
from ..base_classes import DetectorBase, ImagePlaneBase
from ..utils import real_colname, from_currsys, check_keys, interp2


class AutoExposure(Effect):
    """
    Determine DIT and NDIT automatically from ImagePlane

    DIT is determined such that the maximum value in the ``ImagePlane`` fills
    the full well of the detector (``!DET.full_well``) to a given fraction
    (``!OBS.autoexposure.fill_frac``). NDIT is determined such that
    ``DIT`` * ``NDIT`` results in the requested exposure time.

    The requested exposure time is taken from ``!OBS.exptime``.

    The effects sets the parameters `!OBS.dit` and `!OBS.ndit`.

    Example yaml entry
    ------------------
    The parameters `!OBS.exptime`, `!DET.full_well` and `!DET.mindit` should
    be defined as properties in the respective subsections.
    ::
       name: auto_exposure
       description: automatic determination of DIT and NDIT
       class: AutoExposure
       include: True
       kwargs:
           fill_frac: "!OBS.auto_exposure.fill_frac"

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [760]}    # ..todo: change
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ['fill_frac', 'full_well', 'mindit']
        check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, ImagePlaneBase):
            implane_max = np.max(obj.data)
            exptime = from_currsys("!OBS.exptime")
            if exptime is None:
                exptime = from_currsys("!OBS.dit") * from_currsys("!OBS.ndit")
            full_well = from_currsys(self.meta["full_well"])
            fill_frac = from_currsys(self.meta["fill_frac"])
            dit = fill_frac * full_well / implane_max

            # np.ceil so that dit is at most what is required for fill_frac
            ndit = int(np.ceil(exptime / dit))
            dit = exptime / ndit

            # dit must be at least mindit, this might lead to saturation
            # ndit changed so that exptime is not exceeded (hence np.floor)
            if dit < from_currsys(self.meta["mindit"]):
                dit = from_currsys(self.meta["mindit"])
                ndit = int(np.floor(exptime / dit))
                print("Warning: The detector will be saturated!")
                # ..todo: turn into proper warning

            print("Exposure parameters:")
            print("DIT: {:.3f} s     NDIT: {}".format(dit, ndit))

            rc.__currsys__['!OBS.dit'] = dit
            rc.__currsys__['!OBS.ndit'] = ndit

        return obj


class SummedExposure(Effect):
    """
    Simulates a summed stack of ``ndit`` exposures

    """
    def __init__(self, **kwargs):
        super(SummedExposure, self).__init__(**kwargs)
        params = {"z_order": [860]}
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ["dit", "ndit"]
        check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            dit = from_currsys(self.meta["dit"])
            ndit = from_currsys(self.meta["ndit"])

            obj._hdu.data *= dit * ndit

        return obj


class PoorMansHxRGReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(PoorMansHxRGReadoutNoise, self).__init__(**kwargs)
        params = {"z_order": [811],
                  "pedestal_fraction": 0.3,
                  "read_fraction": 0.4,
                  "line_fraction": 0.25,
                  "channel_fraction": 0.05,
                  "random_seed": "!SIM.random.seed",
                  "report_plot_include": False,
                  "report_table_include": False}
        self.meta.update(params)
        self.meta.update(kwargs)

        self.required_keys = ["noise_std", "n_channels", "ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"])
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            from_currsys(self.meta)
            ron_keys = ["noise_std", "n_channels", "channel_fraction",
                        "line_fraction", "pedestal_fraction", "read_fraction"]
            ron_kwargs = {key: self.meta[key] for key in ron_keys}
            ron_kwargs["image_shape"] = det._hdu.data.shape

            ron_frame = make_ron_frame(**ron_kwargs)
            stacked_ron_frame = np.zeros_like(ron_frame)
            for i in range(self.meta["ndit"]):
                dx = np.random.randint(0, ron_frame.shape[1])
                dy = np.random.randint(0, ron_frame.shape[0])
                stacked_ron_frame += np.roll(ron_frame, (dy, dx), axis=(0, 1))

            # .. todo: this .T is ugly. Work out where things are getting switched and remove it!
            det._hdu.data += stacked_ron_frame.T

        return det

    def plot(self, det, new_figure=False, **kwargs):
        import matplotlib.pyplot as plt
        dtcr = self.apply_to(det)
        plt.imshow(dtcr.data, origin="lower")

    def plot_hist(self, det, **kwargs):
        import matplotlib.pyplot as plt
        dtcr = self.apply_to(det)
        plt.hist(dtcr.data.flatten())


class BasicReadoutNoise(Effect):
    def __init__(self, **kwargs):
        super(BasicReadoutNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [811]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

        self.required_keys = ["noise_std", "ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            self.meta = from_currsys(self.meta)
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            noise_std = self.meta["noise_std"] * np.sqrt(float(self.meta["ndit"]))
            det._hdu.data += np.random.normal(loc=0, scale=noise_std,
                                              size=det._hdu.data.shape)

        return det

    def plot(self, det):
        import matplotlib.pyplot as plt
        dtcr = self.apply_to(det)
        plt.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        import matplotlib.pyplot as plt
        dtcr = self.apply_to(det)
        plt.hist(dtcr.data.flatten())


class ShotNoise(Effect):
    def __init__(self, **kwargs):
        super(ShotNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [820]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

    def apply_to(self, det, **kwargs):
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

    def plot(self, det):
        import matplotlib.pyplot as plt
        dtcr = self.apply_to(det)
        plt.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        import matplotlib.pyplot as plt
        dtcr = self.apply_to(det)
        plt.hist(dtcr.data.flatten())


class DarkCurrent(Effect):
    """
    required: dit, ndit, value
    """
    def __init__(self, **kwargs):
        super(DarkCurrent, self).__init__(**kwargs)
        self.meta["z_order"] = [830]

        required_keys = ["value", "dit", "ndit"]
        check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            if isinstance(from_currsys(self.meta["value"]), dict):
                dtcr_id = obj.meta[real_colname("id", obj.meta)]
                dark = from_currsys(self.meta["value"][dtcr_id])
            elif isinstance(from_currsys(self.meta["value"]), float):
                dark = from_currsys(self.meta["value"])
            else:
                raise ValueError("<DarkCurrent>.meta['value'] must be either"
                                 "dict or float: {}".format(self.meta["value"]))

            dit = from_currsys(self.meta["dit"])
            ndit = from_currsys(self.meta["ndit"])

            obj._hdu.data += dark * dit * ndit

        return obj

    def plot(self, det, **kwargs):
        import matplotlib.pyplot as plt
        dit = from_currsys(self.meta["dit"])
        ndit = from_currsys(self.meta["ndit"])
        total_time = dit * ndit
        times = np.linspace(0, 2*total_time, 10)
        dtcr = self.apply_to(det)
        dark_level = dtcr.data[0, 0] / total_time  # just read one pixel
        levels = dark_level * times
        plt.plot(times, levels, **kwargs)
        plt.xlabel("time")
        plt.ylabel("dark level")


class LinearityCurve(Effect):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)
        params = {"z_order": [840],
                  "report_plot_include": True,
                  "report_table_include": False}
        self.meta.update(params)
        self.meta.update(kwargs)

        self.required_keys = ["ndit"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            ndit = from_currsys(self.meta["ndit"])
            incident = self.table["incident"] * ndit
            measured = self.table["measured"] * ndit

            det._hdu.data = np.interp(det._hdu.data, incident, measured)

        return det

    def plot(self, **kwargs):
        import matplotlib.pyplot as plt
        plt.gcf().clf()

        ndit = from_currsys(self.meta["ndit"])
        incident = self.table["incident"] * ndit
        measured = self.table["measured"] * ndit

        plt.loglog(incident, measured, **kwargs)
        plt.xlabel("Incident [ph s$^-1$]")
        plt.ylabel("Measured [e- s$^-1$]")

        return plt.gcf()


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

    def apply_to(self, implane, **kwargs):
        # .. todo: should this be ImagePlaneBase here?
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

    def plot(self, implane, **kwargs):
        import matplotlib.pyplot as plt

        implane = self.apply_to(implane)
        plt.imshow(implane.data, origin="bottom", **kwargs)
        plt.show()


class BinnedImage(Effect):
    def __init__(self, **kwargs):
        super(BinnedImage, self).__init__(**kwargs)
        self.meta["z_order"] = [870]

        self.required_keys = ["bin_size"]
        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            bs = from_currsys(self.meta["bin_size"])
            image = det._hdu.data
            h, w = image.shape
            new_image = image.reshape((h//bs, bs, w//bs, bs))
            det._hdu.data = new_image.sum(axis=3).sum(axis=1)

        return det


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
            dx, dy = min(size[0]-x, n), min(size[1]-y, n)
            image[x:x+dx, y:y+dy] = batch[i:i+dx, j:j+dy]

    return image
