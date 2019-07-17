import numpy as np

from astropy import units as u

from synphot import SourceSpectrum, SpectralElement
from synphot.models import Empirical1D

from .. import rc
from .. import utils
from ..base_classes import SourceBase, ImagePlaneBase
from ..optics.radiometry import RadiometryTable
from .effects import Effect
from .ter_curves import TERCurve


class SurfaceList(Effect):
    """
    A Effect object containing a list of all surfaces in an optical element

    This calls is essentially a wrapper for the RadiometryTable class, which
    contains all the functionality for creating and combining multiple optical
    surfaces


    Parameters
    ----------
    filename : str
        Path to file containing
    table : astropy.Table

    array_dict : dict


    Input format
    ------------

    """

    def __init__(self, **kwargs):
        super(SurfaceList, self).__init__(**kwargs)
        self.meta["z_order"] = []

        self.meta["SIM_MIN_THROUGHPUT"] = rc.__config__["!SIM.spectral.minimum_throughput"]

        self.radiometry_table = RadiometryTable()
        self.radiometry_table.meta.update(self.meta)
        self._emission = None
        self._throughput = None

        data = self.get_data()
        if data is not None:
            self.radiometry_table.add_surface_list(data)

    def apply_to(self, obj, **kwargs):
        if not self.is_empty and isinstance(obj, SourceBase):
            for ii in range(len(obj.spectra)):
                compound_spec = obj.spectra[ii] * self.throughput
                wave = compound_spec.waveset
                spec = compound_spec(wave)
                new_source = SourceSpectrum(Empirical1D, points=wave,
                                            lookup_table=spec)
                obj.spectra[ii] = new_source

        elif not self.is_empty and isinstance(obj, ImagePlaneBase):
            # by calling use_area, the surface area is taken into account, but
            # the units are stuck in PHOTLAM for synphot
            emission = self.get_emission(use_area=True)  # --> PHOTLAM * area
            wave = emission.waveset  # angstrom
            flux = emission(wave)    # PHOTLAM --> ph s-1 cm-2 AA-1 * cm2
            phs = (np.trapz(flux, wave) * u.cm**2).to(u.Unit("ph s-1"))

            obj.hdu.data += phs.value

        return obj

    def fov_grid(self, which="waveset", **kwargs):
        if which == "waveset" and "waverange" in kwargs:
            min_throughput = self.meta["SIM_MIN_THROUGHPUT"]

            wave = np.linspace(min(kwargs["waverange"]),
                               max(kwargs["waverange"]), 100)
            throughput = self.throughput(wave)
            valid_waves = np.where(throughput > min_throughput)[0]
            if len(valid_waves) > 0:
                wave_edges = [min(wave[valid_waves]), max(wave[valid_waves])]
            else:
                raise ValueError("No transmission found above the threshold {} "
                                 "in this wavelength range {}. Did you open "
                                 "the shutter?"
                                 "".format(self.meta["SIM_MIN_THROUGHPUT"],
                                           kwargs["waverange"]))
        else:
            wave_edges = []

        return wave_edges

    def add_surface(self, surface, name, position=-1, add_to_table=True):
        if isinstance(surface, TERCurve):
            surface = surface.surface
        self.radiometry_table.add_surface(surface, name, position, add_to_table)

    def add_surface_list(self, surface_list, prepend=False):
        if isinstance(surface_list, SurfaceList):
            surface_list = surface_list.radiometry_table.table
        elif isinstance(surface_list, RadiometryTable):
            surface_list = surface_list.table

        self.radiometry_table.add_surface_list(surface_list, prepend)

    def get_emission(self, **kwargs):
        if "etendue" in kwargs:
            etendue = kwargs["etendue"]
        elif "etendue" in self.meta:
            etendue = self.meta["etendue"]
        elif "etendue" in self.radiometry_table.meta:
            etendue = self.radiometry_table.meta["etendue"]
        else:
            raise ValueError("etendue must be given in kwargs or .meta")

        return self.radiometry_table.get_emission(etendue=etendue, **kwargs)

    def get_throughput(self, **kwargs):
        return self.radiometry_table.get_throughput(**kwargs)

    def collapse(self, waveset):
        throughput = self.radiometry_table.throughput(waveset)
        self._throughput = SpectralElement(Empirical1D, points=waveset,
                                           lookup_table=throughput)
        emission = self.radiometry_table.emission(waveset)
        self._emission = SourceSpectrum(Empirical1D, points=waveset,
                                        lookup_table=emission)

    @property
    def throughput(self):
        return self.get_throughput()

    @property
    def emission(self):
        return self.get_emission()

    @property
    def area(self):
        if not self.is_empty:
            tbl = self.radiometry_table.table
            outer_col = utils.real_colname("outer", tbl.colnames)
            inner_col = utils.real_colname("inner", tbl.colnames)
            outer = utils.quantity_from_table(outer_col, tbl, u.m)
            inner = utils.quantity_from_table(inner_col, tbl, u.m)
            scope_area = np.max(np.pi / 4 * (outer ** 2 - inner ** 2))
        else:
            scope_area = 0

        return scope_area

    @property
    def is_empty(self):
        return len(self.radiometry_table.table) == 0

