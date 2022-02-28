"""
Influences
Spectral:
    - Red and blue edges of full spectrum
    - Chunks of a large spectral range
Spatial:
    - On-sky borders of Detector Array
    - On-sky borders of Aperture Mask
    - Chunks of large on-sky area
    - Slit mask borders
    - Multiple slits for IFU, mirror-MOS
    - Multiple lenses for fibre-fed MOS

- each Effect should submit a list of volumes

IFU spectroscopy depends on:
    - the tracelist
    - the list of apertures
    - the detector array borders
    - PSF wavelength granularity
    - Atmospheric dispersion

MOS spectroscopy depends on:
    - the tracelist
    - the list of apertures
    - Atmospheric dispersion
    - PSF wavelength granularity
    - the detector array borders

Long slit spectroscopy depends on:
    - the tracelist ra, dec, lam vol --> x, y area
    - the slit aperture
    - the detector array borders
    - PSF wavelength granularity
    - Atmospheric dispersion

Imaging dependent on:
    - Detector array borders
    - PSF wavelength granularity
    - Atmospheric dispersion

"""

from copy import deepcopy, copy
import numpy as np
from astropy.table import Table
from astropy import units as u

from . import fov_manager_utils as fmu
from . import image_plane_utils as ipu
from ..effects import DetectorList
from ..effects import effects_utils as eu
from ..utils import from_currsys

from .fov import FieldOfView
from ..base_classes import FOVSetupBase


class FOVManager:
    """
    A class to manage the (monochromatic) image windows covering the target

    Parameters
    ----------
    effects : list of Effect objects
        Passed from optics_manager.fov_setup_effects

    kwargs
    ------
    All observation parameters as passed from UserCommands

    """
    def __init__(self, effects=[], **kwargs):
        self.meta = {"area": "!TEL.area",
                     "pixel_scale": "!INST.pixel_scale",
                     "plate_scale": "!INST.plate_scale",
                     "wave_min": "!SIM.spectral.wave_min",
                     "wave_mid": "!SIM.spectral.wave_mid",
                     "wave_max": "!SIM.spectral.wave_max",
                     "chunk_size": "!SIM.computing.chunk_size",
                     "max_segment_size": "!SIM.computing.max_segment_size",
                     "sub_pixel": "!SIM.sub_pixel.flag",
                     "sub_pixel_fraction": "!SIM.sub_pixel.fraction",
                     "preload_fovs": "!SIM.computing.preload_field_of_views",
                     "decouple_sky_det_hdrs": "!INST.decouple_detector_from_sky_headers"}
        self.meta.update(kwargs)

        params = from_currsys({"wave_min": self.meta["wave_min"],
                               "wave_max": self.meta["wave_max"]})
        fvl_meta = ["area", "pixel_scale"]
        params["meta"] = from_currsys({key: self.meta[key] for key in fvl_meta})
        self.volumes_list = FovVolumeList(initial_volume=params)

        self.effects = effects
        self._fovs_list = []
        self.is_spectroscope = eu.is_spectroscope(effects)

        if from_currsys(self.meta["preload_fovs"]) is True:
            self._fovs_list = self.generate_fovs_list()

        # ..todo: Make sure this is changed for multiple detectors
        # x_mm, y_mm = ipu.calc_footprint(det_eff.image_plane_header, "D")
        # self.volumes_list.detector_limits = {"xd_min": min(x_mm),
        #                                      "xd_max": max(x_mm),
        #                                      "yd_min": min(y_mm),
        #                                      "yd_max": min(y_mm)}

    def generate_fovs_list(self):
        """
        Generates a series of FieldOfViews objects based self.effects

        Returns
        -------
        fovs : list of FieldOfView objects

        """

        # Ask all the effects to alter the volume_
        params = {"pixel_scale": self.meta["pixel_scale"]}

        for effect in self.effects:
            self.volumes_list = effect.apply_to(self.volumes_list, **params)

        # ..todo: add catch to split volumes larger than chunk_size
        pixel_scale = from_currsys(self.meta["pixel_scale"])
        plate_scale = from_currsys(self.meta["plate_scale"])
        pixel_size = pixel_scale / plate_scale
        plate_scale_deg = plate_scale / 3600.  # ["/mm] / 3600 = [deg/mm]

        chunk_size = from_currsys(self.meta["chunk_size"])
        max_seg_size = from_currsys(self.meta["max_segment_size"])

        split_xs = []
        split_ys = []
        for vol in self.volumes_list:
            vol_pix_area = (vol["x_max"] - vol["x_min"]) * \
                           (vol["y_max"] - vol["y_min"]) / pixel_scale**2
            if vol_pix_area > max_seg_size:
                step = chunk_size * pixel_scale
                split_xs += list(np.arange(vol["x_min"], vol["x_max"], step))
                split_ys += list(np.arange(vol["y_min"], vol["y_max"], step))

        self.volumes_list.split(axis=["x", "y"], value=(split_xs, split_ys))

        fovs = []

        for vol in self.volumes_list:
            xs_min, xs_max = vol["x_min"] / 3600., vol["x_max"] / 3600.
            ys_min, ys_max = vol["y_min"] / 3600., vol["y_max"] / 3600.
            waverange = (vol["wave_min"], vol["wave_max"])
            skyhdr = ipu.header_from_list_of_xy([xs_min, xs_max],
                                                [ys_min, ys_max],
                                                pixel_scale=pixel_scale / 3600.)

            x_sky, y_sky = ipu.calc_footprint(skyhdr)
            x_det = x_sky / plate_scale_deg
            y_det = y_sky / plate_scale_deg
            dethdr = ipu.header_from_list_of_xy(x_det, y_det, pixel_size, "D")
            skyhdr.update(dethdr)

            # useful for spectroscopy mode where slit dimensions is not the same
            # as detector dimensions
            # ..todo: Make sure this changes for multiple image planes
            if from_currsys(self.meta["decouple_sky_det_hdrs"]) is True:
                det_eff = eu.get_all_effects(self.effects, DetectorList)[0]
                dethdr = det_eff.image_plane_header

            fovs += [FieldOfView(skyhdr, waverange, detector_header=dethdr,
                                 **vol["meta"])]

        return fovs

    @property
    def fovs(self):
        if from_currsys(self.meta["preload_fovs"]) is False:
            self._fovs_list = self.generate_fovs_list()
        return self._fovs_list

    @property
    def fov_footprints(self, which="both"):
        return None


class FovVolumeList(FOVSetupBase):
    """
    List of FOV volumes for FOVManager

    Units
    -----
    x, y : [arcsec]
        On-sky coordintates relative to centre of FieldOfView
    wave : [um]
        On-sky wavelengths
    xd, yd : [mm]
        Detector plane coordinated relative to centre of ImagePlane

    """

    def __init__(self, initial_volume={}):

        self.volumes = [{"wave_min": 0.3,
                         "wave_max": 30,
                         "x_min": -1800,
                         "x_max": 1800,
                         "y_min": -1800,
                         "y_max": 1800,
                         "meta": {"area": 0 * u.um**2}
                         }]
        self.volumes[0].update(initial_volume)
        self.detector_limits = {"xd_min": 0,
                                "xd_max": 0,
                                "yd_min": 0,
                                "yd_max": 0}

    def split(self, axis, value):
        """
        Splits the all volumes that include axis=value into two.

        - Loop through all volume dict
        - Find any entries where min < value < max
        - Add two new entries with [min, value], [value, max]

        Parameters
        ----------
        axis : str, list of str
            "wave", "x", "y"
        value : float, list of floats

        Examples
        --------
        ::
            >>> fvl = FovVolumeList()
            >>> fvl.split(axis="wave", value=1.5)
            >>> fvl.split(axis="wave", value=[3.0, 3.1])
            >>> fvl.split(axis=["x", "y"], value=[0, 0])
            >>> fvl.split(axis=["x", "y"], value=([-1, 1], 0))
            >>> fvl.split(axis=["x", "y"], value=([-1, 1], [0, 5]))

        """
        if isinstance(axis, (tuple, list)):
            for ax, val in zip(axis, value):
                self.split(ax, val)
        elif isinstance(value, (tuple, list, np.ndarray)):
            for val in value:
                self.split(axis, val)
        else:
            for i, vol_old in enumerate(self.volumes):
                if vol_old[f"{axis}_min"] < value and vol_old[f"{axis}_max"] > value:
                    vol_new = deepcopy(vol_old)
                    vol_new[f"{axis}_min"] = value
                    vol_old[f"{axis}_max"] = value
                    self.volumes.insert(i+1, vol_new)

    def shrink(self, axis, values):
        """
        - Loop through all volume dict
        - Replace any entries where min < values.min
        - Replace any entries where max > values.max

        Parameters
        ----------
        axis : str
            "wave", "x", "y"
        values : list of 2 floats
            [min, max], [min, None], [None, max]

        Examples
        --------
        ::
            >>> fvl = FovVolumeList()
            >>> fvl.shrink(axis="wave", values=[3.0, 3.1])
            >>> fvl.shrink(axis=["x", "y"], values=([-1, 1], [0, 5]))


        """
        if isinstance(axis, (tuple, list)):
            for ax, val in zip(axis, values):
                self.shrink(ax, val)
        else:
            to_pop = []
            if values[0] is not None:
                for i, vol in enumerate(self.volumes):
                    if vol[f"{axis}_max"] <= values[0]:
                        to_pop += [i]
                    elif vol[f"{axis}_min"] < values[0]:
                        vol[f"{axis}_min"] = values[0]
            if values[1] is not None:
                for i, vol in enumerate(self.volumes):
                    if vol[f"{axis}_min"] >= values[1]:
                        to_pop += [i]
                    if vol[f"{axis}_max"] > values[1]:
                        vol[f"{axis}_max"] = values[1]
            for i in sorted(to_pop)[::-1]:
                self.volumes.pop(i)

    def extract(self, axes, edges):
        """
        Returns new volumes from within all existing volumes

        This method DOES NOT alter the existing self.volumes list
        To include the returned volumes, add them to the self.volumes list

        Parameters
        ----------
        axes : str, list of str
            "wave", "x", "y"
        edges : list, tuple of lists
            Edge points for each axes listed

        Examples
        --------
        ::
            >>> fvl = FovVolumeList()
            >>> fvl.split("x", 0)
            >>> new_vols = fvl.extract(axes=["wave"], edges=([0.5, 0.6]))
            >>> new_vols = fvl.extract(axes=["x", "y"], edges=([-1, 1], [0, 5]))
            >>> new_vols = fvl.extract(axes=["x", "y", "wave"],
            >>>                        edges=([-1, 1], [0, 5], [0.5, 0.6]))
            >>>
            >>> fvl += [new_vols]

        Returns
        -------
        new_vols : list of dicts
            A list of all new volumes extracted from existing volumes

        """
        new_vols = []
        for old_vol in self.volumes:
            add_flag = True
            new_vol = deepcopy(old_vol)
            for axis, edge in zip(axes, edges):
                if edge[0] <= old_vol[f"{axis}_max"] and \
                   edge[1] >= old_vol[f"{axis}_min"]:
                    new_vol[f"{axis}_min"] = max(edge[0], old_vol[f"{axis}_min"])
                    new_vol[f"{axis}_max"] = min(edge[1], old_vol[f"{axis}_max"])
                else:
                    add_flag = False

            if add_flag is True:
                new_vols += [new_vol]

        return new_vols

    def __len__(self):
        return len(self.volumes)

    def __getitem__(self, item):
        return self.volumes[item]

    def __setitem__(self, key, value):
        self.volumes[item] = value

    def __repr__(self):
        text = f"FovVolumeList with [{len(self.volumes)}] volumes:\n"
        for i, vol in enumerate(self.volumes):
            mini_text = ", ".join([f"{key}: {val}" for key, val in vol.items()])
            text += f"  [{i}] {mini_text} \n"

        return text

    def __iadd__(self, other):
        if isinstance(other, list):
            self.volumes += other
        else:
            raise ValueError(f"Can only add lists of volumes: {other}")

        return self

    def __add__(self, other):
        new_self = deepcopy(self)
        new_self += other

        return new_self


# Spectroscopy FOV setup
# ======================
#
# self.effects --> fov_setup_effects [ApertureList, TraceList, DetectorList]
#
# generate a set of sky headers for each aperture in the list
#
# for each trace in tracelist,
#   note the aperture_id
#   note the image_plane_id
#
#   get the pixel_size [from DetectorList.table["pixsize"]]
#   crawl along trace
#       find the red/blue extreme pixel edge wavelengths
#       get the wavelength of the nearest full pixel jump (x or y)
#
#   interpolate the image plane positions of the trace table for the new
#       wavelength array
#   for each new wavelength
#       get the rotation angle of the fov
#       get the shear of the fov wrt the next fov
#
#       generate an image-plane header for each wavelength
#
# get the combined 3d-shifts
# for each trace
#   interpolate the 3d-shifts for all wavelengths
#   apply the required fov-shift to each fov header
