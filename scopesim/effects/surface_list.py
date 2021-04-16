import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt

from .ter_curves import TERCurve
from ..optics import radiometry_utils as rad_utils
from ..optics.surface import PoorMansSurface, SpectralSurface
from ..utils import quantify, from_currsys


class SurfaceList(TERCurve):
    def __init__(self, **kwargs):
        super(SurfaceList, self).__init__(**kwargs)
        params = {"z_order": [20, 120, 520],
                  "minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "etendue": "!TEL.etendue",
                  "report_plot_include": True,
                  "report_table_include": True}
        self.meta.update(params)
        self.meta.update(kwargs)

        self.surfaces = rad_utils.make_surface_dict_from_table(self.table)
        self._surface = None
        self._throughput = None
        self._emission = None

    def fov_grid(self, which="waveset", **kwargs):
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = from_currsys(self.meta)
            wave_min = quantify(self.meta["wave_min"], u.um)
            wave_max = quantify(self.meta["wave_max"], u.um)
            # ..todo:: add 1001 to default.yaml somewhere
            wave = np.linspace(wave_min, wave_max, 1001)
            throughput = self.throughput(wave)
            threshold = self.meta["minimum_throughput"]
            valid_waves = np.where(throughput >= threshold)[0]
            if len(valid_waves) > 0:
                wave_edges = [min(wave[valid_waves]), max(wave[valid_waves])]
            else:
                raise ValueError("No transmission found above the threshold {} "
                                 "in this wavelength range {}. Did you open "
                                 "the shutter?"
                                 "".format(self.meta["minimum_throughput"],
                                           [self.meta["wave_min"],
                                            self.meta["wave_max"]]))
        else:
            wave_edges = []

        return wave_edges

    @property
    def throughput(self):
        if self._throughput is None:
            self._throughput = self.get_throughput()
        return self._throughput

    @property
    def emission(self):
        if "etendue" not in self.meta:
            raise ValueError("self.meta['etendue'] must be set")
        etendue = quantify(from_currsys(self.meta["etendue"]), "m2 arcsec2")
        if self._emission is None:
            self._emission = self.get_emission(etendue)

        return self._emission

    @property
    def surface(self):
        if self._surface is None:
            self._surface = PoorMansSurface(self.emission, self.throughput,
                                            self.meta)

        return self._surface

    @surface.setter
    def surface(self, item):
        self._surface = item

    def get_throughput(self, start=0, end=None, rows=None):
        """ Copied directly from radiometry_table """

        if self.table is None:
            return None
        end = len(self.table) if end is None else end
        end = end + len(self.table) if end < 0 else end
        rows = np.arange(start, end) if rows is None else rows

        thru = rad_utils.combine_throughputs(self.table, self.surfaces, rows)

        return thru

    def get_emission(self, etendue, start=0, end=None, rows=None,
                     use_area=False):
        # ..todo:: work out what this use_area flag means!

        if self.table is None:
            return None
        end = len(self.table) if end is None else end
        end = end + len(self.table) if end < 0 else end
        rows = np.arange(start, end) if rows is None else rows

        emission = rad_utils.combine_emissions(self.table, self.surfaces, rows,
                                               etendue, use_area)

        return emission

    @property
    def area(self):
        areas = [0] + [self.surfaces[key].area.to(u.m**2).value
                       for key in self.surfaces
                       if self.surfaces[key].area is not None]
        return np.max(areas) * u.m**2

    @property
    def is_empty(self):
        return len(self.table) == 0

    # .. todo:: remove: Relic of the old SurfaceList
    def add_surface(self, surface, name=None, position=-1, add_to_table=True):
        if name is None:
            name = surface.meta.get("name", "<unknown surface>")
        if isinstance(surface, TERCurve):
            ter_meta = surface.meta
            surface = surface.surface
            surface.meta.update(ter_meta)

        self.surfaces.update({name: surface})
        self.table = rad_utils.add_surface_to_table(self.table, surface,
                                                    name, position)

    # .. todo:: remove: Relic of the old SurfaceList
    def add_surface_list(self, surface_list, prepend=False):
        if isinstance(surface_list, SurfaceList):
            self.surfaces.update(surface_list.surfaces)
            self.table = rad_utils.combine_tables(surface_list.table,
                                                  self.table, prepend)

    def plot(self, which="x", wavelength=None, ax=None, **kwargs):
        """

        Parameters
        ----------
        which : str
            "x" plots throughput. "t","e","r" plot trans/emission/refl
        wavelength
        kwargs

        Returns
        -------

        """
        import matplotlib.pyplot as plt
        plt.gcf().clf()

        for ii, which_part in enumerate(which):
            ax = plt.subplot(len(which), 1, ii+1)

            # Plot the individual surfaces
            for key in self.surfaces:
                ter = TERCurve(**self.surfaces[key].meta)
                ter.surface = self.surfaces[key]
                ter.plot(which=which_part, wavelength=None, ax=None,
                         new_figure=False,
                         plot_kwargs={"ls": "-", "label": key}, **kwargs)

            # Plot the system surface
            ter = TERCurve(**self.meta)
            ter.surface = self.surface
            ter.plot(which=which_part, wavelength=None, ax=ax, new_figure=False,
                     plot_kwargs={"ls": "-.", "label": "System Throughput"},
                     **kwargs)

        plt.legend()

        return plt.gcf()
