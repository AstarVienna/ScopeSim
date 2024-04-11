# -*- coding: utf-8 -*-
"""
Electronic detector effects - related to detector readout.

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
- UnequalBinnedImage
- Bias - adds constant bias level to readout

Functions:
- make_ron_frame
- pseudo_random_field
"""

import numpy as np

from astropy.io import fits

from .. import rc
from . import Effect
from ..base_classes import DetectorBase, ImagePlaneBase
from ..utils import (from_currsys, figure_factory, check_keys, real_colname,
                     pretty_print_dict, get_logger)


logger = get_logger(__name__)


class DetectorModePropertiesSetter(Effect):
    """
    Set mode specific curr_sys properties for different detector readout modes.

    A little class (``DetectorModePropertiesSetter``) that allows different
    ``"!DET"`` properties to be set on the fly.

    Parameters
    ----------
    mode_properties : dict
        A dictionary containing the DET parameters to be changed for each mode.
        See below for an example yaml entry.

    Examples
    --------
    Add the values for the different detector readout modes to all the relevant
    detector yaml files. In this case the METIS HAWAII (L, M band) and GeoSnap
    (N band) detectors: METIS_DET_IMG_LM.yaml , METIS_DET_IMG_N.yaml
    ::

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
    ::

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

    required_keys = {"mode_properties"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [299, 900]}
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

        self.mode_properties = kwargs["mode_properties"]

    def apply_to(self, obj, **kwargs):
        mode_name = kwargs.get("detector_readout_mode",
                               from_currsys("!OBS.detector_readout_mode",
                                            self.cmds))
        if isinstance(obj, ImagePlaneBase) and mode_name == "auto":
            mode_name = self.select_mode(obj, **kwargs)
            logger.info("Detector mode set to %s", mode_name)

        self.meta["detector_readout_mode"] = mode_name
        props_dict = self.mode_properties[mode_name]
        self.cmds["!OBS.detector_readout_mode"] = mode_name
        for key, value in props_dict.items():
            self.cmds[key] = value

        return obj

    def list_modes(self):
        """Return list of available detector modes."""
        return pretty_print_dict(self.mode_properties)

    def select_mode(self, obj, **kwargs):
        """Automatically select detector mode based on image plane peak value.

        Select the mode with lowest readnoise that does not saturate the
        detector. When all modes saturate, select the mode with the lowest
        saturation level (peak to full_well).
        """
        immax = np.max(obj.data)
        fillfrac = kwargs.get("fill_frac",
                              from_currsys("!OBS.auto_exposure.fill_frac",
                                           self.cmds))

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
    Determine DIT and NDIT automatically from ImagePlane.

    DIT is determined such that the maximum value in the incident photon flux
    (including astronomical source, sky and thermal backgrounds) fills
    the full well of the detector (``!DET.full_well``) to a given fraction
    (``!OBS.autoexposure.fill_frac``). NDIT is determined such that
    ``DIT`` * ``NDIT`` results in the requested exposure time.

    The requested exposure time is taken from ``!OBS.exptime``.

    The effects sets the parameters `!OBS.dit` and `!OBS.ndit`.

    Examples
    --------
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

    required_keys = {"fill_frac", "full_well", "mindit"}

    def __init__(self, **kwargs):
        """
        The effect is the first detector effect, hence essentially operates
        on the `ImagePlane`, mapped to the detector array.
        """
        super().__init__(**kwargs)
        params = {"z_order": [902]}
        self.meta.update(params)
        self.meta.update(kwargs)
        if self.cmds is None:
            logger.error("No cmds present, using default.")
            from scopesim import UserCommands
            self.cmds = UserCommands()

        check_keys(self.meta, self.required_keys, action="error")

    def _dit_above_mindit(self, dit: float) -> bool:
        mindit = from_currsys(self.meta["mindit"], self.cmds)
        if dit < mindit:
            logger.warning("DIT = %.3f s < MINDIT = %.3f s", dit, mindit)
            return False
        return True

    def estimate_dit_ndit(
            self,
            exptime: float,
            image_plane_max: float,
            **kwargs
    ) -> tuple[float, int]:
        """
        Automatically determine DIT and NDIT from exposure time.

        Parameters
        ----------
        exptime : float
            Exposure time in seconds.
        image_plane_max : float
            Maximum pixel value from image plane, used to avoid saturation.

        Returns
        -------
        dit : float
            Detector Integration Time.
        ndit : int
            Number of Integrations.

        """
        # TODO: Remove this silly stuff once currsys works properly...
        full_well = kwargs.get(
            "full_well",
            from_currsys(self.meta["full_well"], self.cmds)
        )
        fill_frac = kwargs.get(
            "fill_frac",
            from_currsys(self.meta["fill_frac"], self.cmds)
        )

        dit_nosat = fill_frac * full_well / image_plane_max
        logger.debug("Required DIT without saturation: %.3f s", dit_nosat)

        # np.ceil so that dit is at most what is required for fill_frac
        ndit = np.ceil(exptime / dit_nosat).astype(int)
        dit = exptime / ndit

        # Note: If the DIT required to avoid saturation is less than MINDIT,
        #       the observation is only possible with likely saturation...
        if not self._dit_above_mindit(dit):
            dit = from_currsys(self.meta["mindit"], self.cmds)
            # NDIT changed so that exptime is not exceeded (hence floor div)
            ndit = max(exptime // dit, 1)
            logger.warning("The detector will likely be saturated!")

        return dit, ndit

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, (ImagePlaneBase, DetectorBase)):
            return obj

        exptime = kwargs.get("exptime",
                             from_currsys("!OBS.exptime", self.cmds))
        mindit = from_currsys(self.meta["mindit"], self.cmds)

        # TODO: Remove this silly try-except once currsys works properly...
        try:
            dit = from_currsys("!OBS.dit", self.cmds)
        except (KeyError, ValueError):
            dit = None
        try:
            ndit = from_currsys("!OBS.ndit", self.cmds)
        except (KeyError, ValueError):
            ndit = None

        if dit and ndit:
            # Both DIT and NDIT are supplied (not None and non-zero), so just
            # use those regardless.
            self.cmds["!OBS.autoexpset"] = False
            # Just log warning in case DIT < MINDIT, don't actually change DIT
            self._dit_above_mindit(dit)
            return obj

        # No DIT or NDIT given, need to determine from exptime
        self.cmds["!OBS.autoexpset"] = True
        if exptime is None:
            logger.warning(
                "Please provide either !OBS.exptime or !OBS.dit + !OBS.ndit")
            if mindit is not None:
                logger.info("Using MINDIT = %.3f s for exposure time.", mindit)
                exptime = mindit
            else:
                logger.warning(
                    "MINDIT not found, falling back to 1 s for exposure time.")
                exptime = 1
        else:
            logger.info("Requested exposure time: %.3f s", exptime)

        dit, ndit = self.estimate_dit_ndit(exptime, obj.data.max(), **kwargs)

        logger.info("Exposure parameters: DIT = %.3f s, NDIT = %d", dit, ndit)
        logger.info("Total exposure time: %.3f s", dit * ndit)

        # TODO: Make sure this goes up far enough in the ChainMap...
        self.cmds["!OBS.dit"] = dit
        self.cmds["!OBS.ndit"] = ndit

        return obj


class SummedExposure(Effect):
    """Simulates a summed stack of ``ndit`` exposures."""

    required_keys = {"dit", "ndit"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [860]}
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            dit = from_currsys(self.meta["dit"], self.cmds)
            ndit = from_currsys(self.meta["ndit"], self.cmds)

            obj._hdu.data *= dit * ndit

        return obj


class Bias(Effect):
    """Adds a constant bias level to readout."""

    required_keys = {"bias"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [855]}
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            biaslevel = from_currsys(self.meta["bias"], self.cmds)
            obj._hdu.data += biaslevel

        return obj

class PoorMansHxRGReadoutNoise(Effect):
    required_keys = {"noise_std", "n_channels", "ndit"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "z_order": [811],
            "pedestal_fraction": 0.3,
            "read_fraction": 0.4,
            "line_fraction": 0.25,
            "channel_fraction": 0.05,
            "random_seed": "!SIM.random.seed",
            "report_plot_include": False,
            "report_table_include": False,
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"],
                                                    self.cmds)
            if self.meta["random_seed"] is not None:
                np.random.seed(self.meta["random_seed"])

            self.meta = from_currsys(self.meta, self.cmds)
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

    def plot(self, det, **kwargs):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data, origin="lower")

    def plot_hist(self, det, **kwargs):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


class BasicReadoutNoise(Effect):
    """Readout noise computed as: ron * sqrt(NDIT)."""

    required_keys = {"noise_std", "ndit"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [811]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            ndit = from_currsys(self.meta["ndit"], self.cmds)
            ron = from_currsys(self.meta["noise_std"], self.cmds)
            noise_std = ron * np.sqrt(float(ndit))

            random_seed = from_currsys(self.meta["random_seed"], self.cmds)
            if random_seed is not None:
                np.random.seed(random_seed)
            det._hdu.data += np.random.normal(loc=0, scale=noise_std,
                                              size=det._hdu.data.shape)

        return det

    def plot(self, det):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


class ShotNoise(Effect):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [820]
        self.meta["random_seed"] = "!SIM.random.seed"
        self.meta.update(kwargs)

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            self.meta["random_seed"] = from_currsys(self.meta["random_seed"],
                                                    self.cmds)
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
                logger.warning(f"Effect ShotNoise: {negvals} negative pixels")
                data[data < 0] = 0

            below = data < 2**20
            above = np.invert(below)
            data[below] = np.random.poisson(data[below]).astype(float)
            data[above] = np.random.normal(data[above], np.sqrt(data[above]))
            new_imagehdu = fits.ImageHDU(data=data, header=det._hdu.header)
            det._hdu = new_imagehdu

        return det

    def plot(self, det):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.imshow(dtcr.data)

    def plot_hist(self, det, **kwargs):
        dtcr = self.apply_to(det)
        fig, ax = figure_factory()
        ax.hist(dtcr.data.flatten())


class DarkCurrent(Effect):
    """
    required: dit, ndit, value
    """

    required_keys = {"value", "dit", "ndit"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [830]

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            if isinstance(from_currsys(self.meta["value"], self.cmds), dict):
                dtcr_id = obj.meta[real_colname("id", obj.meta)]
                dark = from_currsys(self.meta["value"][dtcr_id], self.cmds)
            elif isinstance(from_currsys(self.meta["value"], self.cmds), float):
                dark = from_currsys(self.meta["value"], self.cmds)
            else:
                raise ValueError("<DarkCurrent>.meta['value'] must be either "
                                 f"dict or float, but is {self.meta['value']}")

            dit = from_currsys(self.meta["dit"], self.cmds)
            ndit = from_currsys(self.meta["ndit"], self.cmds)

            obj._hdu.data += dark * dit * ndit

        return obj

    def plot(self, det, **kwargs):
        dit = from_currsys(self.meta["dit"], self.cmds)
        ndit = from_currsys(self.meta["ndit"], self.cmds)
        total_time = dit * ndit
        times = np.linspace(0, 2*total_time, 10)
        dtcr = self.apply_to(det)
        dark_level = dtcr.data[0, 0] / total_time  # just read one pixel
        levels = dark_level * times
        fig, ax = figure_factory()
        ax.plot(times, levels, **kwargs)
        ax.set_xlabel("time")
        ax.set_ylabel("dark level")


class LinearityCurve(Effect):
    """
    Detector linearity effect.

    The detector linearity curve is set in terms of `incident` flux (e/s) and
    `measured` detector values (ADU).

    Examples
    --------

    The effect can be instantiated in various ways.::

        - name: detector_linearity
          class: LinearityCurve
          kwargs:
            filename: FPA_linearity.dat

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

    required_keys = {"ndit"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "z_order": [840],
            "report_plot_include": True,
            "report_table_include": False,
        }
        self.meta.update(params)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            ndit = from_currsys(self.meta["ndit"], self.cmds)
            if self.table is not None:
                incident = self.table["incident"] * ndit
                measured = self.table["measured"] * ndit
            else:
                incident = np.asarray(from_currsys(self.meta["incident"],
                                                   self.cmds)) * ndit
                measured = np.asarray(from_currsys(self.meta["measured"],
                                                   self.cmds)) * ndit
            obj._hdu.data = np.interp(obj._hdu.data, incident, measured)

        return obj

    def plot(self, **kwargs):
        fig, ax = figure_factory()

        ndit = from_currsys(self.meta["ndit"], self.cmds)
        incident = self.table["incident"] * ndit
        measured = self.table["measured"] * ndit

        ax.loglog(incident, measured, **kwargs)
        ax.set_xlabel("Incident [ph s$^-1$]")
        ax.set_ylabel("Measured [e- s$^-1$]")

        return fig


class ReferencePixelBorder(Effect):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
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
        fig, ax = figure_factory()
        ax.imshow(implane.data, origin="bottom", **kwargs)
        # fig.show()


class BinnedImage(Effect):
    required_keys = {"bin_size"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [870]

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            bs = from_currsys(self.meta["bin_size"], self.cmds)
            image = det._hdu.data
            h, w = image.shape
            new_image = image.reshape((h//bs, bs, w//bs, bs))
            det._hdu.data = new_image.sum(axis=3).sum(axis=1)

        return det

class UnequalBinnedImage(Effect):
    required_keys = {"binx","biny"}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta["z_order"] = [870]

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, det, **kwargs):
        if isinstance(det, DetectorBase):
            bx = from_currsys(self.meta["binx"], self.cmds)
            by = from_currsys(self.meta["biny"], self.cmds)
            image = det._hdu.data
            h, w = image.shape
            new_image = image.reshape((h//bx, bx, w//by, by))
            det._hdu.data = new_image.sum(axis=3).sum(axis=1)

        return det


class Quantization(Effect):
    """Converts raw data to whole photons."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "z_order": [825],
            "dtype": "uint32",
        }
        self.meta.update(params)
        self.meta.update(kwargs)

    def _should_apply(self) -> bool:
        if self.cmds is None:
            logger.warning("Cannot access cmds for Quantization effect.")
            return True

        if self.cmds.get("!OBS.autoexpset", False):
            logger.debug("DIT, NDIT determined by AutoExposure -> "
                         "quantization is not applied.")
            return False

        if self.cmds["!OBS.ndit"] > 1:
            logger.debug("NDIT set to 1 -> quantization is not applied.")
            return False

        return True

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, DetectorBase):
            return obj

        if not self._should_apply():
            return obj

        new_dtype = self.meta["dtype"]
        if not np.issubdtype(new_dtype, np.integer):
            logger.warning("Setting quantized data to dtype %s, which is not "
                           "an integer subtype.", new_dtype)

        # This used to create a new ImageHDU with the same header but the data
        # set to the modified data. It should be fine to simply re-assign the
        # data attribute, but just in case it's not...
        logger.debug("Applying quantization to dtype %s.", new_dtype)
        obj._hdu.data = np.floor(obj._hdu.data).astype(new_dtype)

        return obj


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
