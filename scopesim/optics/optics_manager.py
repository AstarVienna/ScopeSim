
from inspect import isclass
from typing import TextIO
from io import StringIO
from collections.abc import Sequence

import numpy as np
from astropy import units as u
from astropy.table import Table
from synphot import SpectralElement, Empirical1D

from .optical_element import OpticalElement
from .. import effects as efs
from ..effects.effects_utils import is_spectroscope
from ..effects.effects_utils import combine_surface_effects
from ..utils import write_report, from_currsys, get_logger
from ..reports.rst_utils import table_to_rst
from .. import rc


logger = get_logger(__name__)


class OpticsManager:
    """
    The workhorse class for dealing with all externally defined Effect objects.

    Parameters
    ----------
    yaml_dicts : list of dict
        The nested dicts describing the Effects from the relevant YAML files,
        which include ``effects`` and ``properties`` sub-dictionaries

    kwargs : expanded dict
        Any extra information not directly related to the optical elements

    """

    def __init__(self, yaml_dicts=None, cmds=None, **kwargs):
        self.optical_elements = []
        self.meta = {}
        self.meta.update(kwargs)
        self._surfaces_table = None
        self._surface_like_effects = None

        self.cmds = cmds
        if self.cmds is None:
            logger.warning("No UserCommands object was passed when initialising OpticsManager")
            self.cmds = rc.__currsys__

        if yaml_dicts is not None:
            self.load_effects(yaml_dicts, **self.meta)

        self.set_derived_parameters()

    def set_derived_parameters(self):

        if "!INST.pixel_scale" not in self.cmds:
            raise ValueError("'!INST.pixel_scale' is missing from the current"
                             "system. Please add this to the instrument (INST)"
                             "properties dict for the system.")
        pixel_scale = self.cmds["!INST.pixel_scale"] << u.arcsec
        if "!TEL.area" in self.cmds and self.cmds["!TEL.area"] != 0:
            area = self.cmds["!TEL.area"] << u.m**2
        else:
            area = self.area
            self.cmds["!TEL.area"] = area
        etendue = area * pixel_scale**2
        self.cmds["!TEL.etendue"] = etendue
        self.cmds["!TEL.area"] = area
        params = {"area": area, "pixel_scale": pixel_scale, "etendue": etendue}
        self.meta.update(params)

    def load_effects(self, yaml_dicts, **kwargs):
        """
        Generate an OpticalElement for each section of the Optical System.

        Make an ``OpticalElement`` for each YAML document in the system.
        For example there should be a YAML document for each of the following:

        - Atmosphere
        - Telescope
        - Relay optics
        - Instrument
        - Detector

        The YAML files can each be separate ``.yaml`` files, or be contained in
        a single ``.yaml`` file separated by a yaml-document-separator:
        ``  \n  ---  \n  ``.

        Parameters
        ----------
        yaml_dicts : list of dicts
            Each YAML dict should contain the descriptions of the Effects
            needed by each ``OpticalElement``.

        """
        if not isinstance(yaml_dicts, Sequence):
            yaml_dicts = [yaml_dicts]
        self.optical_elements.extend(OpticalElement(dic, cmds=self.cmds, **kwargs)
                                     for dic in yaml_dicts if "effects" in dic)

    def add_effect(self, effect, ext=0):
        """
        Add an Effect object to an OpticalElement at index ``ext``.

        Parameters
        ----------
        effect : Effect
            Effect object to be added

        ext : int
            Index number of the desired OpticalElement, contained in the list
            ``self.optical_elements``

        """
        if isinstance(effect, efs.Effect):
            self.optical_elements[ext].add_effect(effect)

    def update(self, **obs_dict):
        """
        Update the meta dictionary with keyword-value pairs.

        Parameters
        ----------
        obs_dict : expanded dict
            Keyword-Value pairs to be added to self.meta

        """
        self.meta.update(**obs_dict)

    def get_all(self, class_type):
        """
        Return list of all effects from all optical elements with `class_type`.

        Parameters
        ----------
        class_type : class object
            The class to be searched for. Must be an class object with
            base-class ``Effect``

        Returns
        -------
        effects : list of Effect objects

        """
        effects = []
        for opt_el in self.optical_elements:
            effects += opt_el.get_all(class_type)

        return effects

    def get_z_order_effects(self, z_level: int):
        """
        Return a list of all effects with a z_order keywords within `z_level`.

        Effect z_order values are classified according to the following:

        - Make a FOV list - z_order = 0..99
        - Make a image plane - z_order = 100..199
        - Apply Source altering effects - z_order = 200..299
        - Apply FOV specific (3D) effects - z_order = 300..399
        - Apply FOV-independent (2D) effects - z_order = 400..499
        - Apply XXX effects - z_order = 500..599
        - Apply XXX effects - z_order = 600..699
        - Apply lambda-independent 2D image plane effects - z_order = 700..799
        - Apply detector effects - z_order = 800..899
        - Apply detector array effects - z_order = 900..999
        - Apply FITS header effects - z_order = 1000...1100

        Parameters
        ----------
        z_level : {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000}
            100-range of z_orders.

        Returns
        -------
        effects : list of Effect objects

        """
        def _gather_effects():
            for opt_el in self.optical_elements:
                yield from opt_el.get_z_order_effects(z_level)

        def _sortkey(eff):
            return next(z % 100 for z in eff.z_order if z >= z_level)

        # return sorted(_gather_effects(), key=_sortkey)
        return list(_gather_effects())

    @property
    def is_spectroscope(self):
        """Return True if any of the effects is a spectroscope."""
        return is_spectroscope(self.all_effects)

    @property
    def image_plane_headers(self):
        """Get headers from detector setup effects."""
        detector_lists = self.detector_setup_effects
        if not detector_lists:
            raise ValueError("No DetectorList objects found.")

        return [det_list.image_plane_header for det_list in detector_lists]

    @property
    def fits_header_effects(self):
        """Get effects with z_order = 1000...1099."""
        return self.get_z_order_effects(1000)

    @property
    def detector_array_effects(self):
        """Get effects with z_order = 900...999."""
        return self.get_z_order_effects(900)

    @property
    def detector_effects(self):
        """Get effects with z_order = 800...899."""
        return self.get_z_order_effects(800)

    @property
    def image_plane_effects(self):
        """Get effects with z_order = 700...799."""
        return self.get_z_order_effects(700)

    @property
    def fov_effects(self):
        """Get effects with z_order = 600...699."""
        return self.get_z_order_effects(600)

    @property
    def source_effects(self):
        """Get effects with z_order = 500...599."""
        return self.get_z_order_effects(500)   # Transmission

    @property
    def detector_setup_effects(self):
        """Get effects with z_order = 400...499 (DetectorLists only!)."""
        return self.get_z_order_effects(400)

    @property
    def image_plane_setup_effects(self):
        """Get effects with z_order = 300...399."""
        return self.get_z_order_effects(300)

    @property
    def fov_setup_effects(self):
        """Get effects with z_order = 200...299."""
        # Working out where to set wave_min, wave_max
        return self.get_z_order_effects(200)

    # TODO: is this ever used anywhere??
    @property
    def surfaces_table(self):
        """Get combined surface table from effects with z_order = 100...199."""
        from copy import deepcopy
        sle_list = self.get_z_order_effects(100)
        sle_list_copy = []
        for eff in sle_list:
            if isinstance(eff, efs.SurfaceList):
                eff_copy = deepcopy(eff)
                eff_copy.table = from_currsys(eff.table, self.cmds)
            else:
                # Avoid infinite recursion in Wheel effects (filter, adc)
                eff_copy = eff
            sle_list_copy.append(eff_copy)

        comb_table = combine_surface_effects(sle_list_copy)
        return comb_table

    @property
    def all_effects(self):
        """Get all effects in all optical elements."""
        return [eff for opt_eff in self.optical_elements for eff in opt_eff]

    @property
    def system_transmission(self):

        wave_unit = u.Unit(from_currsys("!SIM.spectral.wave_unit", self.cmds))
        dwave = from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)
        wave_min = from_currsys("!SIM.spectral.wave_min", self.cmds)
        wave_max = from_currsys("!SIM.spectral.wave_max", self.cmds)
        wave = np.arange(wave_min, wave_max, dwave)
        trans = np.ones_like(wave)
        sys_trans = SpectralElement(Empirical1D, points=wave*u.Unit(wave_unit),
                                    lookup_table=trans)

        for effect in self.get_z_order_effects(100):
            sys_trans *= effect.throughput

        return sys_trans

    @property
    def area(self):
        surf_lists = self.get_all(efs.SurfaceList)
        areas = [0] + [surf_list.area.value for surf_list in surf_lists]
        _area = np.max(areas) * u.m**2

        return _area

    def list_effects(self):
        # unfortunately this does not work because of the astropy internal
        # conversion to np.arrays if all column entries are a list with the
        # same number of entries
        # tbls = [opt_el.list_effects() for opt_el in self.optical_elements]
        # tbl = vstack(tbls)

        # Hence we need to reconstruct the full effects list

        # flat_list = [item for sublist in l for item in sublist]
        all_effs = self.all_effects

        elements = [opt_el.meta["name"] for opt_el in self.optical_elements
                    for eff in opt_el.effects]
        names = [eff.display_name for eff in all_effs]
        classes = [eff.__class__.__name__ for eff in all_effs]
        included = [eff.meta["include"] for eff in all_effs]
        z_orders = [eff.z_order for eff in all_effs]

        colnames = ["element", "name", "class", "included"]     #, "z_orders"
        data = [elements, names, classes, included]             #, z_orders
        data = from_currsys(data, self.cmds)
        tbl = Table(names=colnames, data=data, copy=False)

        return tbl

    def report(self, filename=None, output="rst", rst_title_chars="_^#*+",
               **kwargs):

        rst_str = f"""
List of Optical Elements
{rst_title_chars[0] * 24}

Summary of Effects in Optical Elements:
{rst_title_chars[1] * 39}

.. table::
    :name: tbl:effects_summary

{table_to_rst(self.list_effects(), indent=4)}
"""

        reports = [opt_el.report(rst_title_chars=rst_title_chars[-4:], **kwargs)
                   for opt_el in self.optical_elements]
        rst_str += "\n\n" + "\n\n".join(reports)

        write_report(rst_str, filename, output)

        return rst_str

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        obj = []
        if isclass(item):
            obj += [opt_el.get_all(item) for opt_el in self.optical_elements]
        elif isinstance(item, int):
            obj = self.optical_elements[item]
        elif isinstance(item, str):
            # check for hash-string for getting Effect.meta values
            if item.startswith("#") and "." in item:
                opt_el_name = item.replace("#", "").split(".")[0]
                new_item = item.replace(f"{opt_el_name}.", "")
                obj = self[opt_el_name][new_item]
            else:
                # get all optical elements that match "item"
                obj = [opt_el for opt_el in self.optical_elements
                       if opt_el.meta["name"] == item]

                # add all effects that match "item"
                for opt_el in self.optical_elements:
                    effs = opt_el[item]
                    if not isinstance(effs, list):
                        effs = [effs]
                    obj += effs

        if isinstance(obj, list) and len(obj) == 1:
            obj = obj[0]
        elif isinstance(obj, list) and len(obj) == 0:
            raise ValueError(f"Cannot find object: {item}")

        return obj

    def __setitem__(self, key, value):
        obj = self.__getitem__(key)
        if isinstance(obj, list) and len(obj) > 1:
            logger.warning("%s does not return a singular object:\n %s", key, obj)
        elif isinstance(obj, efs.Effect) and isinstance(value, dict):
            obj.meta.update(value)

    def __contains__(self, key):
        try:
            self[key]
            return True
        except (KeyError, ValueError):
            # FIXME: This should only need KeyError
            return False

    def write_string(self, stream: TextIO) -> None:
        """Write formatted string representation to I/O stream"""
        stream.write(f"{self!s} contains {len(self.optical_elements)} "
                     "OpticalElements\n")
        for opt_elem in enumerate(self.optical_elements):
            opt_elem.write_string(stream, list_effects=False)

    def pretty_str(self) -> str:
        """Return formatted string representation as str."""
        with StringIO() as str_stream:
            self.write_string(str_stream)
            output = str_stream.getvalue()
        return output

    @property
    def display_name(self):
        return self.meta.get("name", self.meta.get("filename", "<empty>"))

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return f"{self.__class__.__name__}: \"{self.display_name}\""
