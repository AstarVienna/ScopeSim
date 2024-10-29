# -*- coding: utf-8 -*-
"""TBA."""

import warnings
from collections import OrderedDict
from copy import deepcopy
from typing import ClassVar

import numpy as np
from astropy import units as u

from .ter_curves import TERCurve, FilterWheelBase
from ..optics import radiometry_utils as rad_utils
from ..optics.surface import PoorMansSurface, SpectralSurface
from ..utils import quantify, from_currsys, figure_factory


class SurfaceList(TERCurve):
    z_order: ClassVar[tuple[int, ...]] = (20, 120, 520)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"minimum_throughput": "!SIM.spectral.minimum_throughput",
                  "etendue": "!TEL.etendue"}
        self.meta.update(params)
        self.meta.update(kwargs)

        self._surface = None
        self._throughput = None
        self._emission = None
        self.surfaces = OrderedDict({})

        if self.table is not None:
            for row in self.table:
                surf_kwargs = deepcopy(self.table.meta)
                surf_kwargs.update(dict(row))
                surf_kwargs["cmds"] = self.cmds
                surf_kwargs["filename"] = from_currsys(surf_kwargs["filename"], self.cmds)
                self.surfaces[surf_kwargs["name"]] = SpectralSurface(**surf_kwargs)

    def fov_grid(self, which="waveset", **kwargs):
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release.", DeprecationWarning, stacklevel=2)
        wave_edges = []
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = from_currsys(self.meta, self.cmds)
            wave_min = quantify(self.meta["wave_min"], u.um)
            wave_max = quantify(self.meta["wave_max"], u.um)
            # ..todo:: add 1001 to default.yaml somewhere
            wave = np.linspace(wave_min, wave_max, 1001)
            throughput = self.throughput(wave)
            threshold = self.meta["minimum_throughput"]
            valid_waves = np.where(throughput >= threshold)[0]

            if not len(valid_waves):
                msg = ("No transmission found above the threshold "
                       f"{self.meta['minimum_throughput']} in this wavelength "
                       f"range {[self.meta['wave_min'], self.meta['wave_max']]}."
                       " Did you open the shutter?")
                raise ValueError(msg)

            wave_edges = [min(wave[valid_waves]), max(wave[valid_waves])]
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
        etendue = quantify(from_currsys(self.meta["etendue"], self.cmds),
                           "m2 arcsec2")
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
        """Copied directly from radiometry_table."""
        if self.table is None:
            return None
        if end is None:
            end = len(self.table)
        if end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return rad_utils.combine_throughputs(self.table, self.surfaces, rows)

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
        if isinstance(surface, (TERCurve, FilterWheelBase)):
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
            # new_tbl = from_currsys(surface_list.table, self.cmds),
            self.table = rad_utils.combine_tables(surface_list.table, self.table, prepend)

    def plot(self, which="x", wavelength=None, *, axes=None, **kwargs):
        """Plot TER curves.

        Parameters
        ----------
        which : {"x", "t", "e", "r"}, optional
            "x" plots throughput. "t","e","r" plot trans/emission/refl.
            Can be a combination, e.g. "tr" or "tex" to plot each.
        wavelength : array_like, optional
            Passed to TERCurve.plot() for each surface. The default is None.
        axes : matplotlib axes, optional
            If given, plot into existing axes. The default is None.

        Returns
        -------
        fig : matplotlib figure
            Figure containing plots.

        """
        if axes is None:
            fig, axes = figure_factory(len(which), 1, iterable_axes=True)
        else:
            fig = axes.figure
            self._axes_guard(which, axes)

        for ter, ax in zip(which, axes):
            # Plot the individual surfaces
            # TODO: do we want separate plots for these? if yes, how (row/col)?
            for key, surface in self.surfaces.items():
                curve = TERCurve(**surface.meta)
                curve.surface = surface
                kwargs.update(plot_kwargs={"ls": "-", "label": key})
                curve.plot(ter, wavelength, axes=ax, **kwargs)

            # Plot the system surface
            # TODO: self is a subclass of TERCurve, why create again??
            curve = TERCurve(**self.meta)
            curve.surface = self.surface
            kwargs.update(plot_kwargs={"ls": "-.",
                                       "label": "System Throughput"})
            curve.plot(ter, axes=ax, **kwargs)

        fig.legend()

        return fig
