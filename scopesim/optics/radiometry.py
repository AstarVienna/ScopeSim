from collections import OrderedDict

import numpy as np
from astropy import units as u

from ..utils import real_colname, quantify
from .radiometry_utils import combine_emissions, combine_throughputs, \
    combine_tables, add_surface_to_table, add_surface_to_dict, \
    make_surface_from_row


class RadiometryTable:
    def __init__(self, tables=(), **kwargs):
        self.meta = {"area": None}
        self.meta.update(kwargs)

        self.table = None
        self.surfaces = OrderedDict({})

        if len(tables) > 0:
            self.add_surface_list(tables)

    def add_surface_list(self, surface_list, prepend=False):
        self.table = combine_tables(surface_list, self.table, prepend=prepend)

        r_name = real_colname("name", self.table.colnames)
        for row in self.table:
            if row[r_name] not in self.surfaces:
                surf = make_surface_from_row(row)
                self.add_surface(surf, row[r_name], position=-1,
                                 add_to_table=False)

    def add_surface(self, surface, name, position=-1, add_to_table=True):
        if self.table is None:
            raise ValueError("Cannot add surface without <self>.table template."
                             "Please add an empty table to define column names")

        if position < 0:
            position += len(self.table) + 1
        self.surfaces = add_surface_to_dict(self.surfaces, surface,
                                            name, position)
        if add_to_table:
            self.table = add_surface_to_table(self.table, surface,
                                              name, position)

    def get_throughput(self, start=0, end=None, rows=None):

        if self.table is None:
            return None

        if end is None:
            end = len(self.table)
        elif end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return combine_throughputs(self.table, self.surfaces, rows)

    def get_emission(self, etendue, start=0, end=None, rows=None,
                     use_area=False):
        if self.table is None:
            return None

        if end is None:
            end = len(self.table)
        elif end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return combine_emissions(self.table, self.surfaces, rows, etendue,
                                 use_area)

    @property
    def emission(self):
        if "etendue" not in self.meta:
            raise ValueError("self.meta['etendue'] must be set")
        etendue = quantify(self.meta["etendue"], "m2 arcsec2")

        return self.get_emission(etendue)

    @property
    def throughput(self):
        return self.get_throughput()

    def plot(self, what="all", rows=None):
        raise NotImplemented

    def __getitem__(self, item):
        return self.surfaces[item]

    def __repr__(self):
        return self.table.__repr__()


