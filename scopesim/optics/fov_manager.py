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

from copy import deepcopy
from typing import TextIO
from io import StringIO
from collections.abc import Iterable, Iterator, MutableSequence
from itertools import chain

import numpy as np
from astropy import units as u
from astropy.wcs import WCS

from . import image_plane_utils as ipu
from ..effects import DetectorList
from ..effects import effects_utils as eu
from ..utils import from_currsys, get_logger

from .fov import FieldOfView
from ..base_classes import FOVSetupBase


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
            # ..todo: Make sure this changes for multiple image planes
            if from_currsys(self.meta["decouple_sky_det_hdrs"], self.cmds):
                det_eff = eu.get_all_effects(self.effects, DetectorList)[0]
                dethdr = det_eff.image_plane_header

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


class FovVolumeList(FOVSetupBase, MutableSequence):
    """
    List of FOV volumes for FOVManager.

    Units
    -----
    x, y : [arcsec]
        On-sky coordintates relative to centre of FieldOfView
    wave : [um]
        On-sky wavelengths
    xd, yd : [mm]
        Detector plane coordinates relative to centre of ImagePlane

    """

    def __init__(self, initial_volume=None):
        if initial_volume is None:
            initial_volume = {}

        self.volumes = [{
            "wave_min": 0.3,
            "wave_max": 30,
            # 100 square degree should be enough for everyone!
            # Survey telescopes can have a large Field of View though:
            # - OmegaCAM / KiDS: 1 sq. degree
            # - DES: 4 sq. degree
            # - Pan-STARRS 7 sq. degree
            # - DREAMS: 6 sq. degree
            # TODO: Why not put the entire sky here?
            "x_min": -1800 * 10,
            "x_max": 1800 * 10,
            "y_min": -1800 * 10,
            "y_max": 1800 * 10,
            "meta": {
                "area": 0 * u.um**2,
                "aperture_id": 0,
            },
        }]
        self.volumes[0].update(initial_volume)  # .. TODO: Careful! This overwrites meta
        self.detector_limits = {
            "xd_min": 0,
            "xd_max": 0,
            "yd_min": 0,
            "yd_max": 0,
        }

    def split(self, axis, value, aperture_id=None) -> None:
        """
        Split the all volumes that include axis=value into two.

        - Loop through all volume dict
        - Find any entries where min < value < max
        - Add two new entries with [min, value], [value, max]

        Parameters
        ----------
        axis : {"wave", "x", "y"}, or list thereof
            Which axis (``str``) or axes (``list[str]``) to use.
        value : float, list of floats
        aperture_id : int, optional
            Default None. If ``None``, split all volumes. If ``int``, only split
            volumes with this ``aperture_id`` in the meta dict

        Examples
        --------
        ::
            >>> fvl = FovVolumeList()
            >>> fvl.split(axis="wave", value=1.5)
            >>> fvl.split(axis="wave", value=[3.0, 3.1])
            >>> fvl.split(axis=["x", "y"], value=[0, 0])
            >>> fvl.split(axis=["x", "y"], value=([-1, 1], 0))
            >>> fvl.split(axis=["x", "y"], value=([-1, 1], [0, 5]))
            >>> fvl.split(axis="wave", value=3.0, aperture_id=1)
            >>> fvl.split(axis="wave", value=3.0, aperture_id=None)

        """
        if isinstance(axis, Iterable) and not isinstance(axis, str):
            for ax, val in zip(axis, value):
                self.split(ax, val)
            return

        if isinstance(value, Iterable):
            for val in value:
                self.split(axis, val)
            return

        for vol in self:
            if _chk_ap_id(aperture_id, vol):
                continue
            if vol[f"{axis}_min"] >= value or vol[f"{axis}_max"] <= value:
                continue
            new_vol = deepcopy(vol)
            new_vol[f"{axis}_min"] = value
            vol[f"{axis}_max"] = value
            self.insert(self.index(vol) + 1, new_vol)

    def shrink(self, axis, values, aperture_id=None) -> None:
        """
        Trim axes to new min/max value(s).

        - Loop through all volume dict
        - Replace any entries where min < values.min
        - Replace any entries where max > values.max

        Parameters
        ----------
        axis : {"wave", "x", "y"} or list thereof
            Which axis (``str``) or axes (``list[str]``) to use.
        values : list of 2 floats
            [min, max], [min, None], [None, max]
        aperture_id : int, optional
            Default None. If ``None``, shrink all volumes. If ``int``, only
            shrink volumes with this ``aperture_id`` in the meta dict

        Examples
        --------
        ::
            >>> fvl = FovVolumeList()
            >>> fvl.shrink(axis="wave", values=[3.0, 3.1])
            >>> fvl.shrink(axis="wave", values=[2.9, 3.1], aperture_id=1)
            >>> fvl.shrink(axis=["x", "y"], values=([-1, 1], [0, 5]))


        """
        # FIXME: Isn't this method just the same as setting self.volumes to the
        #        output list of self.extract()?? Except the None values.
        if isinstance(axis, (tuple, list)):
            for ax, val in zip(axis, values):
                self.shrink(ax, val)
            return

        to_pop = []

        for vol in self:
            if _chk_ap_id(aperture_id, vol):
                continue

            if values[0] is not None:
                if vol[f"{axis}_max"] <= values[0]:
                    to_pop.append(self.index(vol))
                    continue
                vol[f"{axis}_min"] = max(values[0], vol[f"{axis}_min"])

            if values[1] is not None:
                if vol[f"{axis}_min"] >= values[1]:
                    to_pop.append(self.index(vol))
                    continue
                vol[f"{axis}_max"] = min(values[1], vol[f"{axis}_max"])

        for idx in reversed(sorted(to_pop)):
            self.pop(idx)

    def extract(self, axes, edges, aperture_id=None):
        """
        Return new volumes from within all existing volumes.

        This method DOES NOT alter the existing self.volumes list
        To include the returned volumes, add them to the self.volumes list

        Parameters
        ----------
        axes : list of either {"wave", "x", "y"}
            Which axis (``list`` of single ``str``) or axes (``list[str]``)
            to use. Must be ``list`` in either case.
        edges : list, tuple of lists
            Edge points for each axes listed
        aperture_id : int, optional
            Default None. If ``None``, extract from all volumes. If ``int``,
            only extract from volumes with this `aperture_id` in the meta dict

        Examples
        --------
        ::
            >>> fvl = FovVolumeList()
            >>> fvl.split("x", 0)
            >>> new_vols = fvl.extract(axes=["wave"], edges=([0.5, 0.6], ))
            >>> new_vols = fvl.extract(axes=["x", "y"], edges=([-1, 1], [0, 5]))
            >>> new_vols = fvl.extract(axes=["x", "y", "wave"],
            >>>                        edges=([-1, 1], [0, 5], [0.5, 0.6]))
            >>> new_vols = fvl.extract(axes=["x", "y"], edges=([-1, 1], [0, 5]),
            >>>                        aperture_id=1)
            >>>
            >>> fvl += [new_vols]

        Returns
        -------
        new_vols : list of dicts
            A list of all new volumes extracted from existing volumes

        """
        def _get_new_vols():
            for vol in self:
                if _chk_ap_id(aperture_id, vol):
                    continue
                if not all(_volume_in_range(vol, axis, edge) for axis, edge
                           in zip(axes, edges)):
                    continue
                new_vol = deepcopy(vol)
                for axis, edge in zip(axes, edges):
                    new_vol[f"{axis}_min"] = max(edge[0], vol[f"{axis}_min"])
                    new_vol[f"{axis}_max"] = min(edge[1], vol[f"{axis}_max"])
                yield new_vol

        return list(_get_new_vols())

    def __len__(self):
        return len(self.volumes)

    def __getitem__(self, index):
        return self.volumes[index]

    def __setitem__(self, index, value):
        self.volumes[index] = value

    def __delitem__(self, index):
        del self.volumes[index]

    def insert(self, index, value):
        self.volumes.insert(index, value)

    def write_string(self, stream: TextIO) -> None:
        """Write formatted string representation to I/O stream."""
        n_vol = len(self.volumes)
        stream.write(f"FovVolumeList with {n_vol} volume{'s'*(n_vol>1)}:")
        max_digits = len(str(n_vol))

        for i_vol, vol in enumerate(self.volumes):
            pre = "\n└─" if i_vol == n_vol - 1 else "\n├─"
            stream.write(f"{pre}[{i_vol:>{max_digits}}]:")

            pre = "\n  " if i_vol == n_vol - 1 else "\n│ "
            n_key = len(vol)
            for i_key, (key, val) in enumerate(vol.items()):
                subpre = "└─" if i_key == n_key - 1 else "├─"
                stream.write(f"{pre}{subpre}{key}: {val}")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.volumes[0]})"

    def __str__(self) -> str:
        with StringIO() as str_stream:
            self.write_string(str_stream)
            output = str_stream.getvalue()
        return output

    def __add__(self, other):
        # TODO: Is this functionality actually used anywhere?
        new_self = deepcopy(self)
        new_self += other

        return new_self

    def _repr_pretty_(self, p, cycle):
        """For ipython."""
        if cycle:
            p.text(f"{self.__class__.__name__}(...)")
        else:
            p.text(str(self))


def _volume_in_range(vol: dict, axis: str, edge) -> bool:
    return edge[0] <= vol[f"{axis}_max"] and edge[1] >= vol[f"{axis}_min"]


def _chk_ap_id(ap_id, vol) -> bool:
    return ap_id is not None and ap_id != vol["meta"]["aperture_id"]


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
