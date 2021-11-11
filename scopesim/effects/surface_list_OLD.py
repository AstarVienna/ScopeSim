
import numpy as np

from astropy import units as u

from synphot import SourceSpectrum, SpectralElement
from synphot.models import Empirical1D

from .. import rc
from .. import utils
from ..utils import quantify
from ..base_classes import SourceBase, ImagePlaneBase, FieldOfViewBase
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
        self.meta["z_order"] = [20, 120]
        self.meta["minimum_throughput"] = "!SIM.spectral.minimum_throughput"
        self.meta["wave_min"] = "!SIM.spectral.wave_min"
        self.meta["wave_max"] = "!SIM.spectral.wave_max"
        self.meta.update(kwargs)

        self.radiometry_table = RadiometryTable()
        self.radiometry_table.meta.update(self.meta)
        self._emission = None
        self._throughput = None

        data = self.get_data()
        if data is not None:
            self.radiometry_table.add_surface_list(data)

    def apply_to(self, obj, **kwargs):
        """
        obj == SourceBase - applies throughput
        obj == ImagePlaneBase - applies emission if Imager
        obj == FieldOfViewBase - applies emission if Spectrograph

        """
        if isinstance(obj, SourceBase) and not self.is_empty:
            self.meta = utils.from_currsys(self.meta)
            for ii in range(len(obj.spectra)):
                spec = obj.spectra[ii]
                wave_val = spec.waveset.value
                wave_unit = spec.waveset.unit  # angstrom
                wave_min = quantify(self.meta["wave_min"], u.um).to(u.AA)
                wave_max = quantify(self.meta["wave_max"], u.um).to(u.AA)
                mask = (wave_val > wave_min.value) * (wave_val < wave_max.value)

                wave = ([wave_min.value] +
                        list(wave_val[mask]) +
                        [wave_max.value]) * wave_unit
                thru = self.throughput(wave)
                flux = spec(wave)
                flux *= thru
                new_source = SourceSpectrum(Empirical1D, points=wave,
                                            lookup_table=flux)
                obj.spectra[ii] = new_source

        elif isinstance(obj, ImagePlaneBase) and not self.is_empty:
            # by calling use_area, the surface area is taken into account, but
            # the units are stuck in PHOTLAM for synphot
            emission = self.get_emission(use_area=True)  # --> PHOTLAM * area
            if emission is not None:
                wave = emission.waveset  # angstrom
                flux = emission(wave)    # PHOTLAM --> ph s-1 cm-2 AA-1 * cm2
                phs = (np.trapz(flux, wave) * u.cm**2).to(u.Unit("ph s-1"))
            else:
                phs = 0 * (u.ph / u.s)

            obj.hdu.data += phs.value

        elif isinstance(obj, FieldOfViewBase) and not self.is_empty:
            # ..todo:: Super hacky, FIX THIS!!
            emission = self.get_emission(use_area=True)  # --> PHOTLAM * area
            if emission is not None:
                wave_val = emission.waveset.value
                wave_unit = emission.waveset.unit   # angstrom
                wave_min = quantify(obj.meta["wave_min"], u.um).to(wave_unit)
                wave_max = quantify(obj.meta["wave_max"], u.um).to(wave_unit)
                mask = (wave_val > wave_min.value) * (wave_val < wave_max.value)

                wave = ([wave_min.value] + list(wave_val[mask]) +
                        [wave_max.value]) * wave_unit
                flux = emission(wave)    # PHOTLAM --> ph s-1 cm-2 AA-1 * cm2
                phs = (np.trapz(flux, wave) * u.cm**2).to(u.Unit("ph s-1"))
            else:
                phs = 0 * (u.ph / u.s)

            obj.hdu.data += phs.value

        return obj

    def fov_grid(self, which="waveset", **kwargs):
        if which == "waveset":
            self.meta.update(kwargs)
            self.meta = utils.from_currsys(self.meta)
            wave_min = utils.quantify(self.meta["wave_min"], u.um)
            wave_max = utils.quantify(self.meta["wave_max"], u.um)
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

    def add_surface(self, surface, name, position=-1, add_to_table=True):
        if isinstance(surface, TERCurve):
            ter_meta = surface.meta
            surface = surface.surface
            surface.meta.update(ter_meta)
        self.radiometry_table.add_surface(surface, name, position, add_to_table)

    def add_surface_list(self, surface_list, prepend=False):
        if isinstance(surface_list, SurfaceList):
            surface_list = surface_list.radiometry_table.table
        elif isinstance(surface_list, RadiometryTable):
            surface_list = surface_list.table

        self.radiometry_table.add_surface_list(surface_list, prepend)

    def get_emission(self, **kwargs):
        etendue = rc.__currsys__["!TEL.etendue"]
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
