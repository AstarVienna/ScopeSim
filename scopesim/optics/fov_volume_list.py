# -*- coding: utf-8 -*-
"""Outsourced from `fov_manager` to avoid circular imports."""


from copy import deepcopy
from typing import TextIO
from io import StringIO
from collections.abc import Iterable, MutableSequence

from astropy import units as u

from ..utils import get_logger


logger = get_logger(__name__)


class FovVolumeList(MutableSequence):
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
