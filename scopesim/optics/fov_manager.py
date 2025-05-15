# -*- coding: utf-8 -*-
"""
Influences.

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

from collections.abc import Iterator
from itertools import chain

import numpy as np
from astropy.wcs import WCS

from . import image_plane_utils as ipu
from ..effects import DetectorList
from ..effects import effects_utils as eu
from ..utils import from_currsys, get_logger

from .fov import FieldOfView
from .fov_volume_list import FovVolumeList


logger = get_logger(__name__)


class FOVManager:
    """
    A class to manage the (monochromatic) image windows covering the target.

    Parameters
    ----------
    effects : list of Effect objects
        Passed from optics_manager.fov_setup_effects

    kwargs
    ------
    All observation parameters as passed from UserCommands

    """

    def __init__(self, effects=None, cmds=None, **kwargs):
        self.meta = {
            "area": "!TEL.area",
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
            "decouple_sky_det_hdrs": "!INST.decouple_detector_from_sky_headers",
            "aperture_id": 0,
        }
        self.meta.update(kwargs)
        self.cmds = cmds

        params = from_currsys({"wave_min": self.meta["wave_min"],
                               "wave_max": self.meta["wave_max"]},
                              self.cmds)
        fvl_meta = ["area", "pixel_scale", "aperture_id"]
        params["meta"] = from_currsys({key: self.meta[key] for key in fvl_meta},
                                      self.cmds)
        self.volumes_list = FovVolumeList(initial_volume=params)

        self.effects = effects or []
        self._fovs_list = []
        self.is_spectroscope = eu.is_spectroscope(self.effects)

        if from_currsys(self.meta["preload_fovs"], self.cmds):
            logger.debug("Generating initial fovs_list.")
            self._fovs_list = list(self.generate_fovs_list())

    def _get_splits(self, pixel_scale):
        chunk_size = from_currsys(self.meta["chunk_size"], self.cmds)
        max_seg_size = from_currsys(self.meta["max_segment_size"], self.cmds)

        for vol in self.volumes_list:
            vol_pix_area = ((vol["x_max"] - vol["x_min"]) *
                            (vol["y_max"] - vol["y_min"]) / pixel_scale**2)
            if vol_pix_area > max_seg_size:
                step = chunk_size * pixel_scale

                # These are not always integers, unlike in the tests.
                # See for example HAWKI/test_hawki/test_full_package_hawki.py.
                # The np.arange can therefore not be changed to just a range.
                yield (np.arange(vol["x_min"], vol["x_max"], step),
                       np.arange(vol["y_min"], vol["y_max"], step))

    def generate_fovs_list(self) -> Iterator[FieldOfView]:
        """
        Generate a series of FieldOfViews objects based self.effects.

        Yields
        ------
        Iterator[FieldOfView]
            Generator-Iterator of FieldOfView objects.

        """
        # TODO: The generator is currently always converted to a list, but that
        #       might not be necessary. Investigate in the future...

        # Ask all the effects to alter the volume_
        params = {"pixel_scale": self.meta["pixel_scale"]}

        for effect in self.effects:
            self.volumes_list = effect.apply_to(self.volumes_list, **params)

        # ..todo: add catch to split volumes larger than chunk_size
        pixel_scale = from_currsys(self.meta["pixel_scale"], self.cmds)
        plate_scale = from_currsys(self.meta["plate_scale"], self.cmds)

        splits = (chain.from_iterable(split)
                  for split in zip(*self._get_splits(pixel_scale)))

        self.volumes_list.split(axis=["x", "y"], value=splits)

        for vol in self.volumes_list:
            xs_min, xs_max = vol["x_min"] / 3600., vol["x_max"] / 3600.
            ys_min, ys_max = vol["y_min"] / 3600., vol["y_max"] / 3600.
            waverange = (vol["wave_min"], vol["wave_max"])
            skyhdr = ipu.header_from_list_of_xy([xs_min, xs_max],
                                                [ys_min, ys_max],
                                                pixel_scale=pixel_scale / 3600.)

            dethdr, _ = ipu.det_wcs_from_sky_wcs(
                WCS(skyhdr), pixel_scale, plate_scale)
            skyhdr.update(dethdr.to_header())

            # useful for spectroscopy mode where slit dimensions is not the same
            # as detector dimensions
            # TODO: Make sure this changes for multiple image planes
            if from_currsys(self.meta["decouple_sky_det_hdrs"], self.cmds):
                det_eff = eu.get_all_effects(self.effects, DetectorList)[0]
                dethdr = det_eff.image_plane_header
                # TODO: Why is this .image_plane_header and not
                #       .detector_headers()[0] or something?

            yield FieldOfView(skyhdr,
                              waverange,
                              detector_header=dethdr,
                              cmds=self.cmds,
                              **vol["meta"])

    @property
    def fovs(self):
        # There two lines were not here before #258, but somehow not including
        # them will mess things up as FOVs multipy like rabbits...
        # Should investigate why at some point...
        if self._fovs_list:
            logger.debug("Returning existing fovs_list.")
            return self._fovs_list

        if not from_currsys(self.meta["preload_fovs"], self.cmds):
            logger.debug("Generating new fovs_list.")
            self._fovs_list = list(self.generate_fovs_list())
        return self._fovs_list

    @property
    def fov_footprints(self):
        return None


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
