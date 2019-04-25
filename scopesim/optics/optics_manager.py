import warnings

from astropy import units as u

from .. import effects as efs
from ..effects.effects_utils import combine_surface_effects
from .optical_element import OpticalElement


class OpticsManager:
    """
    The workhorse class for dealing with all externally defined Effect objects

    Parameters
    ----------
    yaml_dicts : list of dict
        The nested dicts describing the Effects from the relevant YAML files

    kwargs : **dict
        Any keyword-value pairs from a config file


    """
    def __init__(self, yaml_dicts=[], **kwargs):
        self.optical_elements = [OpticalElement({"name": "misc"})]
        self.meta = {}
        self.meta.update(kwargs)
        self._surfaces_table = None

        if yaml_dicts is not None:
            self.load_effects(yaml_dicts)

        self.set_pixel_scale(yaml_dicts)

    def set_pixel_scale(self, yaml_dicts):
        # .. todo: hack. Think up a more elegant way to store SIM_PIXEL_SCALE

        if isinstance(yaml_dicts, dict):
            yaml_dicts = [yaml_dicts]

        if "SIM_PIXEL_SCALE" not in self.meta:
            self.meta["SIM_PIXEL_SCALE"] = None

        if self.meta["SIM_PIXEL_SCALE"] is None:
            pixel_scale = None
            for yaml_dict in yaml_dicts:
                if "properties" in yaml_dict:
                    if "pixel_scale" in yaml_dict["properties"]:
                        pixel_scale = yaml_dict["properties"]["pixel_scale"]

            self.meta["SIM_PIXEL_SCALE"] = pixel_scale

    def load_effects(self, yaml_dicts):
        """
        Generate an OpticalElement for each section of the Optical System

        Make an OpticalElement for each YAML document in the system. For example
        there should be a YAML document for each of the following:

        - Atmosphere
        - Telescope
        - Relay optics
        - Instrument
        - Detector

        The YAML files can each be separate .yaml files, or be contained in a
        single .yaml file separated by a yaml-document-separator: ``\n --- \n``

        Parameters
        ----------
        yaml_dicts : list of dicts
            Each YAML dict should contain the descriptions of the Effects needed
            by each OpticalElement

        """

        if isinstance(yaml_dicts, dict):
            yaml_dicts = [yaml_dicts]
        self.optical_elements += [OpticalElement(dic) for dic in yaml_dicts]

    def add_effect(self, effect, ext=0):
        """
        Add an Effect object to an OpticalElement at index ``ext``

        Parameters
        ----------
        effect : Effect
            Effect object to be added

        ext : int
            Index number of the desired OpticalElement, contained in the list
            self.optical_elements

        """
        if isinstance(effect, efs.Effect):
            self.optical_elements[ext].add_effect(effect)

    def update(self, **obs_dict):
        """
        Update the meta dictionary with keyword-value pairs

        Parameters
        ----------
        obs_dict : **dict
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
        z_level : int
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
    def image_plane_header(self):
        detector_lists = self.get_all(efs.DetectorList)
        header = detector_lists[0].image_plane_header

        if len(detector_lists) != 1:
            warnings.warn("None or more than one DetectorList found. Using the"
                          " first instance.{}".format(detector_lists))

        return header

    @property
    def detector_effects(self):
        dtcr_effects = []
        for opt_el in self.optical_elements:
            dtcr_effects += opt_el.get_z_order_effects([500, 599])

        return dtcr_effects

    @property
    def image_plane_effects(self):
        imp_effects = [self.surfaces_table]
        for opt_el in self.optical_elements:
            imp_effects += opt_el.get_z_order_effects([400, 499])

        return imp_effects

    @property
    def fov_effects(self):
        fov_effects = []
        for opt_el in self.optical_elements:
            fov_effects += opt_el.get_z_order_effects([300, 399])

        return fov_effects

    @property
    def source_effects(self):
        src_effects = [self.surfaces_table]
        for opt_el in self.optical_elements:
            src_effects += opt_el.get_z_order_effects([200, 299])

        return src_effects

    @property
    def image_plane_setup_effects(self):
        implane_setup_effects = []
        for opt_el in self.optical_elements:
            implane_setup_effects += opt_el.get_z_order_effects([100, 199])

        return implane_setup_effects

    @property
    def fov_setup_effects(self):
        fovmanager_effects = [self.surfaces_table]
        for opt_el in self.optical_elements:
            fovmanager_effects += opt_el.get_z_order_effects([0, 99])

        return fovmanager_effects

    @property
    def surfaces_table(self):
        surface_like_effects = []
        for opt_el in self.optical_elements:
            surface_like_effects += opt_el.ter_list

        surf_table = None
        if len(surface_like_effects) > 0:
            pixel_scale = self.meta["SIM_PIXEL_SCALE"] * u.arcsec
            surf_table = combine_surface_effects(surface_like_effects)
            surf_table.meta["etendue"] = surf_table.area * pixel_scale**2

            self._surfaces_table = surf_table

        return surf_table

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        if isinstance(item, efs.Effect):
            effects = []
            for opt_el in self.optical_elements:
                effects += opt_el.get_all(item)
            return effects
        elif isinstance(item, int):
            return self.optical_elements[item]

    def __repr__(self):
        msg = "\nOpticsManager contains {} OpticalElements \n" \
              "".format(len(self.optical_elements))
        for ii, opt_el in enumerate(self.optical_elements):
            msg += '[{}] "{}" contains {} effects \n' \
                   ''.format(ii, opt_el.meta["name"], len(opt_el.effects))

        return msg
