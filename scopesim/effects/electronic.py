"""
Electronic detector effects - related to detector readout

Classes:
- DetectorModePropertiesSetter - set parameters for readout mode
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
import logging

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from .. import rc
from . import Effect
from ..base_classes import DetectorBase, ImagePlaneBase
from ..utils import from_currsys
from .. import utils


class DetectorModePropertiesSetter(Effect):
    """
    Sets mode specific curr_sys properties for different detector readout modes

    A little class (DetectorModePropertiesSetter) that allows different "!DET"
    properties to be set on the fly.

    Parameters
    ----------
    mode_properties : dict
        A dictionary containing the DET parameters to be changed for each mode
        See below for an example yaml entry.

    Examples
    --------

    Add the values for the different detector readout modes to all the relevant
    detector yaml files. In this case the METIS HAWAII (L, M band) and GeoSnap
    (N band) detectors: METIS_DET_IMG_LM.yaml , METIS_DET_IMG_N.yaml

        - name: lm_detector_readout_parameters
          class: DetectorModePropertiesSetter
          kwargs:
            mode_properties:
              fast:
                mindit: 0.04
                full_well: !!float 1e5
                ron: 70
              slow:
                mindit: 1.3
                full_well: !!float 1e5
                ron: 14

    Add the OBS dict entry !OBS.detector_readout_mode to the properties section
    of the mode_yamls descriptions in the default.yaml files.

        mode_yamls:
          - object: observation
            alias: OBS
            name: lss_l
            yamls:
              ...
            properties:
              ...
              detector_readout_mode: slow

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [299, 900]}
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ['mode_properties']
        utils.check_keys(self.meta, required_keys, action="error")

        self.mode_properties = kwargs['mode_properties']

    def apply_to(self, obj, **kwargs):
        mode_name = kwargs.get('detector_readout_mode',
                               from_currsys("!OBS.detector_readout_mode"))
        if isinstance(obj, ImagePlaneBase) and mode_name == "auto":
            mode_name = self.select_mode(obj, **kwargs)
            print("Detector mode set to", mode_name)

        self.meta['detector_readout_mode'] = mode_name
        props_dict = self.mode_properties[mode_name]
        rc.__currsys__["!OBS.detector_readout_mode"] = mode_name
        for key, value in props_dict.items():
            rc.__currsys__[key] = value

        return obj

    def list_modes(self):
        """Return list of available detector modes"""
        return utils.pretty_print_dict(self.mode_properties)

    def select_mode(self, obj, **kwargs):
        """Automatically select detector mode based on image plane peak value

        Select the mode with lowest readnoise that does not saturate the detector.
        When all modes saturate, select the mode with the lowest saturation level
        (peak to full_well).
        """
        immax = np.max(obj.data)
        fillfrac = kwargs.get("fill_frac",
                              from_currsys("!OBS.auto_exposure.fill_frac"))

        goodmodes = []
        goodron = []
        badmodes = []
        satlevel = []
        for modeid, modeprops in self.mode_properties.items():
            mindit = modeprops["!DET.mindit"]
            fullwell = modeprops["!DET.full_well"]
            if immax * mindit < fillfrac * fullwell:
                goodmodes.append(modeid)
                goodron.append(modeprops["!DET.readout_noise"])
            else:
                badmodes.append(modeid)
                satlevel.append(immax * mindit / fullwell)

        if not goodmodes:
            return badmodes[np.argmin(satlevel)]

        return goodmodes[np.argmin(goodron)]


class AutoExposure(Effect):
    """
    Determine DIT and NDIT automatically from ImagePlane

    DIT is determined such that the maximum value in the incident photon flux
    (including astronomical source, sky and thermal backgrounds) fills
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
        """
        The effect is the first detector effect, hence essentially operates
        on the `ImagePlane`, mapped to the detector array.
        """
        super().__init__(**kwargs)
        params = {"z_order": [902]}
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ['fill_frac', 'full_well', 'mindit']
        utils.check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, (ImagePlaneBase, DetectorBase)):
            implane_max = np.max(obj.data)
            exptime = kwargs.get('exptime', from_currsys("!OBS.exptime"))
            mindit = from_currsys(self.meta["mindit"])

            if exptime is None:
                exptime = from_currsys("!OBS.dit") * from_currsys("!OBS.ndit")
            print(f"Requested exposure time: {exptime:.3f} s")

            if exptime < mindit:
                print(f"    increased to MINDIT: {mindit:.3f} s")
                exptime = mindit

            full_well = from_currsys(self.meta["full_well"])
            fill_frac = kwargs.get("fill_frac",
                                   from_currsys(self.meta["fill_frac"]))
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
            print(f"                DIT: {dit:.3f} s  NDIT: {ndit}")
            print(f"Total exposure time: {dit * ndit:.3f} s")

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
        utils.check_keys(self.meta, required_keys, action="error")

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
        utils.check_keys(self.meta, self.required_keys, action="error")

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
        dtcr = self.apply_to(det)
        plt.imshow(dtcr.data, origin="lower")

    def plot_hist(self, det, **kwargs):
        dtcr = self.apply_to(det)
        plt.hist(dtcr.data.flatten())


class BasicReadoutNoise(Effect):
    """Readout noise computed as: ron * sqrt(NDIT)"""
    def __init__(self, **kwargs):
        super(BasicReadoutNoise, self).__init__(**kwargs)
        self.meta["z_order"] = [811]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

        self.required_keys = ["noise_std", "ndit"]
        utils.check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            ndit = from_currsys(self.meta["ndit"])
            ron = from_currsys(self.meta["noise_std"])
            noise_std = ron * np.sqrt(float(ndit))

            random_seed = from_currsys(self.meta["random_seed"])
            if random_seed is not None:
                np.random.seed(random_seed)
            det._hdu.data += np.random.normal(loc=0, scale=noise_std,
                                              size=det._hdu.data.shape)

        return det

    def plot(self, det):
        dtcr = self.apply_to(det)
        plt.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
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
            # Windows has a problem with generating poisson values above 2**30
            # Above ~100 counts the poisson and normal distribution are
            # basically the same. For large arrays the normal distribution
            # takes only 60% as long as the poisson distribution
            data = det._hdu.data

            # Check if there are negative values in the data
            negvals = np.sum(data < 0)
            if negvals:
                logging.warning(f"Effect ShotNoise: {negvals} negative pixels")

            below = data < 2**20
            above = np.invert(below)
            data[below] = np.random.poisson(data[below]).astype(float)
            data[above] = np.random.normal(data[above], np.sqrt(data[above]))
            data = np.floor(data)
            new_imagehdu = fits.ImageHDU(data=data, header=det._hdu.header)
            det._hdu = new_imagehdu

        return det

    def plot(self, det):
        dtcr = self.apply_to(det)
        plt.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
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
        utils.check_keys(self.meta, required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            if isinstance(from_currsys(self.meta["value"]), dict):
                dtcr_id = obj.meta[utils.real_colname("id", obj.meta)]
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
    """
    Detector linearity effect

    The detector linearity curve is set in terms of `incident` flux (e/s) and `measured`
    detector values (ADU).

    Examples
    --------

    The effect can be instantiated in various ways.
    - name: detector_linearity
      class: LinearityCurve
      kwargs:
        file_name: FPA_linearity.dat

    - name: detector_linearity
      class: LinearityCurve
      kwargs:
        array_dict: {incident: [0, 77000, 999999999999],
                     measured: [0, 77000, 77000]}

    - name: detector_linearity
      class: LinearityCurve
      kwargs:
        incident: [0, 77000, 99999999]
        measured: [0, 77000, 77000]
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [840],
                  "report_plot_include": True,
                  "report_table_include": False}
        self.meta.update(params)
        self.meta.update(kwargs)

        self.required_keys = ["ndit"]
        utils.check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            ndit = from_currsys(self.meta["ndit"])
            if self.table is not None:
                incident = self.table["incident"] * ndit
                measured = self.table["measured"] * ndit
            else:
                incident = np.asarray(from_currsys(self.meta["incident"])) * ndit
                measured = np.asarray(from_currsys(self.meta["measured"])) * ndit
            obj._hdu.data = np.interp(obj._hdu.data, incident, measured)

        return obj

    def plot(self, **kwargs):
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
        implane = self.apply_to(implane)
        plt.imshow(implane.data, origin="bottom", **kwargs)
        plt.show()


class BinnedImage(Effect):
    def __init__(self, **kwargs):
        super(BinnedImage, self).__init__(**kwargs)
        self.meta["z_order"] = [870]

        self.required_keys = ["bin_size"]
        utils.check_keys(self.meta, self.required_keys, action="error")

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
