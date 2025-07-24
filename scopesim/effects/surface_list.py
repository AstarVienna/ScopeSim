# -*- coding: utf-8 -*-
"""TBA."""

import warnings
from collections import OrderedDict
from copy import deepcopy
from typing import ClassVar

import numpy as np
from astropy import units as u

from .ter_curves import TERCurve
from ..optics.surface import PoorMansSurface, SpectralSurface
from ..utils import quantify, from_currsys, figure_factory, real_colname


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

        return self.combine_throughputs(rows)

    def get_emission(self, etendue, start=0, end=None, rows=None,
                     use_area=False):
        # TODO: deprecate use_area flag

        if self.table is None:
            return None
        end = len(self.table) if end is None else end
        end = end + len(self.table) if end < 0 else end
        rows = np.arange(start, end) if rows is None else rows

        emission = self.combine_emissions(rows)

        return emission

    def combine_throughputs(self, rows_indexes):
        if len(self.table) == 0:
            return None

        r_name = real_colname("name", self.table.colnames)
        r_action = real_colname("action", self.table.colnames)

        throughput = None
        for ii, row_num in enumerate(rows_indexes):

            row = self.table[row_num]
            surf = self.surfaces[row[r_name]]
            action_attr = row[r_action]
            if action_attr == "":
                raise ValueError(f"No action in surf.meta: {surf.meta}")

            if isinstance(surf, SpectralSurface):
                surf_throughput = getattr(surf, action_attr)

                if ii == 0:
                    throughput = deepcopy(surf_throughput)
                else:
                    throughput = throughput * surf_throughput

        return throughput

    def combine_emissions(self, row_indexes):
        """
        Combine thermal emission from a series of surfaces.

        The function traces thermal emission through an optical system, taking
        into account the finite reflectivities/transmissivities and emissivities
        of the surfaces. The function assumes that etendue is conserved through
        the system, i.e. surfaces are neither over- nor undersized.

        Parameters
        ----------
        tbl : astropy Table
            Required columns are `name` and `action` (reflection or transmission)
        surfaces: OrderedDict of SpectralSurface
            Keys are the names from tbl, values are of type `SpectralSurface`
        row_indexes : list of int
            Rows of tbl (i.e. surfaces) to combine

        Returns
        -------
        SourceSpectrum

        """
        if len(self.table) == 0:
            return None

        r_name = real_colname("name", self.table.colnames)
        r_action = real_colname("action", self.table.colnames)

        emission = None
        for row_num in row_indexes:
            row = self.table[row_num]
            surf = self.surfaces[row[r_name]]
            action_attr = row[r_action]

            if isinstance(surf, SpectralSurface):
                surf_throughput = getattr(surf, action_attr)
                if emission is not None:
                    emission = emission * surf_throughput

                if emission is None:
                    emission = deepcopy(surf.emission)
                else:
                    emission = emission + surf.emission

        return emission

    @property
    def area(self):
        areas = [0] + [self.surfaces[key].area.to_value(u.m**2)
                       for key in self.surfaces
                       if self.surfaces[key].area is not None]
        return np.max(areas) * u.m**2

    @property
    def is_empty(self):
        return len(self.table) == 0

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
