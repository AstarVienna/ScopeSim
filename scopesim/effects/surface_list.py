import numpy as np

from .effects import Effect
from ..optics import radiometry_utils as rad_utils
from ..utils import quantify, from_currsys


class SurfaceList(Effect):
    def __init__(self, **kwargs):
        super(SurfaceList, self).__init__(**kwargs)
        params = {"z_order": [20, 120, 520],
                  "minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "wave_min": "!SIM.spectral.wave_min",
                  "wave_max": "!SIM.spectral.wave_max",
                  "wave_bin": "!SIM.spectral.spectral_resolution",
                  "wave_unit": "!SIM.spectral.wave_unit",
                  "etendue": "!TEL.etendue",
                  "area": "!TEL.area"}
        self.meta.update(params)
        self.meta.update(kwargs)

        self.surfaces = rad_utils.make_surface_dict_from_table(self.table)
        self._surface = None
        self._throughput = None

    @property
    def throughput(self):
        if self._throughput is None:
            self._throughput = self.get_throughput()
        return self._throughput

    def get_throughput(self, start=0, end=None, rows=None):
        """ Copied directly from radiometry_table """

        if self.table is None:
            return None

        if end is None:
            end = len(self.table)
        elif end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return rad_utils.combine_throughputs(self.table, self.surfaces, rows)

    @property
    def emission(self):
        if "etendue" not in self.meta:
            raise ValueError("self.meta['etendue'] must be set")
        etendue = quantify(from_currsys(self.meta["etendue"]), "m2 arcsec2")

        return self.get_emission(etendue)

    def get_emission(self, etendue, start=0, end=None, rows=None,
                     use_area=False):
        # ..todo:: work out what this use_area flag means!

        if self.table is None:
            return None

        if end is None:
            end = len(self.table)
        elif end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return rad_utils.combine_emissions(self.table, self.surfaces, rows,
                                           etendue, use_area)
