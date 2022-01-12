import logging
from inspect import isclass

import numpy as np
from astropy import units as u
from astropy.table import Table
from synphot import SpectralElement, Empirical1D

from .optical_element import OpticalElement
from .. import effects as efs
from ..effects.effects_utils import is_spectroscope
from ..effects.effects_utils import combine_surface_effects
from ..utils import write_report, from_currsys
from ..reports.rst_utils import table_to_rst
from .. import rc


class OpticsManager:
    """
    The workhorse class for dealing with all externally defined Effect objects

    Parameters
    ----------
    yaml_dicts : list of dict
        The nested dicts describing the Effects from the relevant YAML files,
        which include ``effects`` and ``properties`` sub-dictionaries

    kwargs : expanded dict
        Any extra information not directly related to the optical elements

    """

    def __init__(self, yaml_dicts=[], **kwargs):
        self.optical_elements = []
        self.meta = {}
        self.meta.update(kwargs)
        self._surfaces_table = None
        self._surface_like_effects = None

        if yaml_dicts is not None:
            self.load_effects(yaml_dicts, **self.meta)

        self.set_derived_parameters()

    def set_derived_parameters(self):

        if "!INST.pixel_scale" not in rc.__currsys__:
            raise ValueError("!INST.pixel_scale is missing from the current"
                             "system. Please add this to the instrument (INST)"
                             "properties dict for the system.")
        pixel_scale = rc.__currsys__["!INST.pixel_scale"] * u.arcsec
        area = self.area
        etendue = area * pixel_scale**2
        rc.__currsys__["!TEL.etendue"] = etendue
        rc.__currsys__["!TEL.area"] = area

        params = {"area": area, "pixel_scale": pixel_scale, "etendue": etendue}
        self.meta.update(params)

    def load_effects(self, yaml_dicts, **kwargs):
        """
        Generate an OpticalElement for each section of the Optical System

        Make an ``OpticalElement`` for each YAML document in the system.
        For example there should be a YAML document for each of the following:

        - Atmosphere
        - Telescope
        - Relay optics
        - Instrument
        - Detector

        The YAML files can each be separate ``.yaml`` files, or be contained in
        a single ``.yaml`` file separated by a yaml-document-separator:
        ``\n --- \n``

        Parameters
        ----------
        yaml_dicts : list of dicts
            Each YAML dict should contain the descriptions of the Effects needed
            by each ``OpticalElement``

        """

        if isinstance(yaml_dicts, dict):
            yaml_dicts = [yaml_dicts]
        self.optical_elements += [OpticalElement(dic, **kwargs)
                                  for dic in yaml_dicts if "effects" in dic]

    def add_effect(self, effect, ext=0):
        """
        Add an Effect object to an OpticalElement at index ``ext``

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
        Update the meta dictionary with keyword-value pairs

        Parameters
        ----------
        obs_dict : expanded dict
            Keyword-Value pairs to be added to self.meta

        """
        self.meta.update(**obs_dict)

    def get_all(self, class_type):
        """
        Return a list of all effects from all optical elements with `class_type`

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

    def get_z_order_effects(self, z_level):
        """
        Return a list of all effects with a z_order keywords within z_level

        Effect z_order values are classified according to the following:

        - Make a FOV list - z_order = 0..99
        - Make a image plane - z_order = 100..199
        - Apply Source altering effects - z_order = 200..299
        - Apply FOV specific (3D) effects - z_order = 300..399
        - Apply FOV-independent (2D) effects - z_order = 400..499

        Parameters
        ----------
        z_level : int, tuple
            [0, 100, 200, 300, 400, 500]

        Returns
        -------
        effects : list of Effect objects

        """

        effects = []
        for opt_el in self.optical_elements:
            effects += opt_el.get_z_order_effects(z_level)

        return effects

    @property
    def is_spectroscope(self):
        return bool(len(self.get_all(efs.SpectralTraceList)))

    @property
    def image_plane_headers(self):
        detector_lists = self.detector_setup_effects
        headers = [det_list.image_plane_header for det_list in detector_lists]

        if len(detector_lists) == 0:
            raise ValueError("No DetectorList objects found. {}"
                             "".format(detector_lists))

        return headers

    @property
    def detector_array_effects(self):
        return self.get_z_order_effects(900)

    @property
    def detector_effects(self):
        return self.get_z_order_effects(800)

    @property
    def image_plane_effects(self):
        effects = self.get_z_order_effects(700)
        return effects

    @property
    def fov_effects(self):
        effects = self.get_z_order_effects(600)
        return effects

    @property
    def source_effects(self):
        return self.get_z_order_effects(500)   # Transmission

    @property
    def detector_setup_effects(self):
        # !!! Only DetectorLists go in here !!!
        return self.get_z_order_effects(400)

    @property
    def image_plane_setup_effects(self):
        return self.get_z_order_effects(300)

    @property
    def fov_setup_effects(self):
        # Working out where to set wave_min, wave_max
        return self.get_z_order_effects(200)

    @property
    def surfaces_table(self):
        if self._surfaces_table is None:
            surface_like_effects = self.get_z_order_effects(100)
            self._surfaces_table = combine_surface_effects(surface_like_effects)
        return self._surfaces_table

    @property
    def system_transmission(self):

        wave_unit = u.Unit(rc.__currsys__["!SIM.spectral.wave_unit"])
        dwave = rc.__currsys__["!SIM.spectral.spectral_bin_width"]
        wave_min = rc.__currsys__["!SIM.spectral.wave_min"]
        wave_max = rc.__currsys__["!SIM.spectral.wave_max"]
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
        all_effs = [eff for opt_eff in self.optical_elements for eff in opt_eff]

        elements = [opt_el.meta["name"] for opt_el in self.optical_elements
                    for eff in opt_el.effects]
        names = [eff.display_name for eff in all_effs]
        classes = [eff.__class__.__name__ for eff in all_effs]
        included = [eff.meta["include"] for eff in all_effs]
        z_orders = [eff.meta["z_order"] for eff in all_effs]

        colnames = ["element", "name", "class", "included"]     #, "z_orders"
        data = [elements, names, classes, included]             #, z_orders
        data = from_currsys(data)
        tbl = Table(names=colnames, data=data, copy=False)

        return tbl

    def report(self, filename=None, output="rst", rst_title_chars="_^#*+",
               **kwargs):

        rst_str = """
List of Optical Elements
{}

Summary of Effects in Optical Elements:
{}

.. table::
    :name: tbl:effects_summary

{}
""".format(rst_title_chars[0] * 24,
           rst_title_chars[1] * 39,
           table_to_rst(self.list_effects(), indent=4))

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
            obj = [opt_el for opt_el in self.optical_elements
                   if opt_el.meta["name"] == item]
            for opt_el in self.optical_elements:
                obj += opt_el[item]

        if isinstance(obj, list) and len(obj) == 1:
            obj = obj[0]
        elif isinstance(obj, list) and len(obj) == 0:
            raise ValueError(f"Cannot find object: {item}")

        return obj

    def __setitem__(self, key, value):
        obj = self.__getitem__(key)
        if isinstance(obj, list) and len(obj) > 1:
            logging.warning("{} does not return a singular object:\n {}"
                          "".format(key, obj))
        elif isinstance(obj, efs.Effect) and isinstance(value, dict):
            obj.meta.update(value)

    def __repr__(self):
        msg = "\nOpticsManager contains {} OpticalElements \n" \
              "".format(len(self.optical_elements))
        for ii, opt_el in enumerate(self.optical_elements):
            msg += '[{}] "{}" contains {} effects \n' \
                   ''.format(ii, opt_el.meta["name"], len(opt_el.effects))

        return msg

    def __str__(self):
        name = self.meta.get("name", self.meta.get("filename", "<empty>"))
        return '{}: "{}"'.format(type(self).__name__, name)
