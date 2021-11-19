# """
# This module contains classes which control how a simulation is run
#
# Module Summary
# --------------
# UserCommands is essentially a dictionary that holds all the variables that
# the user may wish to change. It also has some set variables like ``pix_res``
# that can be accessed directly, instead of from the dictionary.
#
# UserCommands is imported directly into the simmetis package and is accessible
# from the main package - ``simmetis.UserCommands``
#
# Classes
# -------
# ``UserCommands(filename, sim_data_dir=<path to data files>)``
#
# Routines
# --------
#
# * ``dump_defaults(filename="./", selection="freq")``
# * ``dump_chip_layout(dir="./")``
#
#
# See Also
# --------
# Classes that require a ``UserCommands`` object directly include:
#
# * ``Detector``
# * ``OpticalTrain``
#
# Examples
# --------
# By default ``UserCommands`` contains the parameters needed to generate the MICADO
# optical train:
#
#     >>> my_cmds = scopesim.UserCommands()
#
# To list the keywords that are available:
#
#     >>> my_cmds.keys()
#     ...
#
#
# The UserCommands object also contains smaller dictionaries for each category of
# keywords - e.g. for the keywords for the instrument:
#
#     >>> my_cmds.inst
#     ...
#
# """

import os
import glob
import logging
import copy

import yaml
import numpy as np
import astropy.io.ascii as ioascii

from scopesim import rc

from scopesim.utils import find_file, zendist2airmass, \
                    airmass2zendist
from scopesim.effects import atmospheric_refraction
from scopesim import server as svr

from OLD_code import OLD_user_commands_utils as cutils, OLD_spectral as sc

__all__ = ["UserCommands"]


class UserCommands(object):
    """
    An extended dictionary with the parameters needed for running a simulation

    Summary
    -------
    A :class:`.UserCommands` object contains a dictionary which holds all the
    keywords from the ``default.config`` file. It also has attributes which
    represent the frequently used variables, i.e. ``pix_res``, ``lam_bin_edges``
    , ``exptime``, etc

    ``<UserCommands>.cmds`` is a dictionary that holds all the variables
    the user may wish to change. It also has some set variables like
    ``<UserCommands>.pix_res`` that can be accessed directly, instead of from
    the dictionary.

    ``UserCommands`` is imported directly into the scopesim package and is
    accessable from the main package - ``scopesim.UserCommands``

    If UserCommands is called without any arguments, the default values for
    MICADO and the E-ELT are used.


    Parameters
    ----------
    filename : str, optional
        path to the user's .config file

    Attributes
    ----------
    Internal dictionaries
    cmds : dict (collections.OrderedDict)
        the dictionary which holds all the keyword-value pairs needed for
        running a simualtion
    obs : dict (collections.OrderedDict)
        parameters about the observation
    sim : dict (collections.OrderedDict)
        parameters about the simualtion
    atmo : dict (collections.OrderedDict)
        parameters about the atmosphere
    scope : dict (collections.OrderedDict)
        parameters about the telescope
    inst : dict (collections.OrderedDict)
        parameters about the instrument
    fpa : dict (collections.OrderedDict)
        parameters about the detector array (FPA - Focal Plane Array)
    hxrg : dict (collections.OrderedDict)
        parameters about the chip noise (HxRG - HAWAII 4RG chip series)

    Attributes pertaining to the purely spectral data sets (e.g. transmission
    curves, stellar spectra, etc)
    lam : np.ndarray
        a vector containing the centres of the wavelength bins used when
        resampling the spectra or transmission curves
    lam_res : float
        [um] the resolution of the ``lam``

    Attributes pertaining to the binning in spectral space for which different
    PSFs need to be used
    lam_psf_res : float
        [um] the spectal "distance" between layers - i.e. width of the bins
    lam_bin_edges : array-like
        [um] the edge of the spectral bin for each layer
    lam_bin_centers : array-like
        [um] the centres of the spectral bin for each layer

    Attributes pertaining to the binning in the spatial plane
    pix_res : float
        [arcsec] the internal (oversampled) spatial resolution of the simulation
    fpa_res : float
        [arcsec] the field of view of the individual pixels on the detector

    General attributes
    verbose : bool
        Flag for printing intermediate results to the screen (default=True)
    exptime : float
        [s] exposure time of a single DIT
    diameter : float
        [m] outer diamter of the primary aperture (i.e. M1)
    area : float
        [m^2] effective area of the primary aperture (i.e. M1)
    filter : str
        [BVRIzYJHKKs,user] filter used for the observation

    Methods
    -------
    update(new_dict)
        updates the current ``UserCommands`` object from another dict-like object
    keys()
        returns the keys in the ``UserCommands.cmds`` dictionary
    values()
        returns the values in the ``UserCommands.cmds`` dictionary


    Examples
    --------
    By default ``UserCommands`` contains the parameters needed to generate the
    MICADO optical train::

        >>> import scopesim
        >>> empty_cmds = scopesim.UserCommands()
        >>> micado_cmds = scopesim.UserCommands(instrument="MICADO")
        >>> mcao_wide_cmds = scopesim.UserCommands(instrument="MICADO",
        ...                                        mode="MODE_MCAO_WIDE")

    To list the keywords that are available and their values::

        >>> mcao_wide_cmds.keys
        >>> mcao_wide_cmds.values

    The ``UserCommands`` object also contains smaller dictionaries for each
    category of keywords - e.g. for the keywords describing the instrument::

        >>> mcao_wide_cmds.atmo
        >>> mcao_wide_cmds.inst

    """

    def __init__(self, filename=None, sim_data_dir=None,
                 instrument=None, mode=None, filter_name=None):

        """
        Create an extended dictionary of simulation parameters

        Parameters
        ----------
        filename : str, optional
            path to the user's .config file
        sim_data_dir : str, optional
            path to directory where instrument data are stored

        SimCADO needs to know where the instrument-specific data files
        are stored. This can be specified in the config file (keyword
        `SIM_DATA_DIR`) or in the call to `UserCommands`. The default
        configuration file does not include a valid sim_data_dir, hence at
        least one of the parameters `filename` or `sim_data_dir` must be
        provided.

        """

        self.pkg_dir = rc.__pkg_dir__
        self.data_dir = rc.__data_dir__

        # read in the default keywords
        # Don't use self.update because we need to add all the valid keywords
        self.cmds = copy.deepcopy(rc.__old_config__)
        self.cmds.update(rc.__rc__)
        self._convert_python_types()

        self.cmds["CONFIG_USER"] = filename
        self.cmds["CONFIG_DEFAULT"] = os.path.join(rc.__pkg_dir__,
                                                   "OLD.default.config")
        # read in the users wishes
        # Use self.update so that we reject all the invalid keywords
        if filename is not None:
            self.update(filename)

        # option sim_data_dir overrides values in config files
        if sim_data_dir is not None:
            self.cmds['SIM_DATA_DIR'] = sim_data_dir

        if self.cmds["SIM_OVERSAMPLING"] == "yes":
            self.cmds["PSF_MODE"] = "oversample"
        else:
            self.cmds["PSF_MODE"] = "linear_interp"

        if instrument is not None:
            self.set_instrument(instrument)

        if mode is not None:
            self.set_mode(mode)

        if filter_name is not None:
            self.select_filter(filter_name)

        if self.verbose and filename is not None:
            print("Read in parameters from " + filename)

    def update(self, new_input):
        """
        Update multiple entries of a ``UserCommands`` dictionary

        ``update(new_dict)`` takes either a normal python ``dict`` object or a
        ``UserCommands`` object. Only keywords that match those in the
        ``UserCommands`` object will be updated. Misspelled keywords raise an
        error.

        To update single items in the dictionary, it is recommended to simply
        call the key and update the value - i.e ``<UserCommands>[key] = value``.

        Parameters
        ----------
        new_input : str, dict, ``UserCommands``, file path


        Raises
        ------
        KeyError
            If a parameter is not found in ``self.cmds``.

        Examples
        --------
        View the default tests_commands::

            >>> import scopesim
            >>> my_cmds = scopesim.UserCommands()
            >>> print(my_cmds.cmds)

        Change a single command::

            >>> my_cmds["OBS_EXPTIME"] = 60

        Change a series of tests_commands at once::

            >>> new_cmds = {"OBS_EXPTIME" : 60 , "OBS_NDIT" : 10}
            >>> my_cmds.update(new_cmds)

        .. todo::
            Add the new functionality to docs:
            - can update with strings, dicts, filenames, or UserCommands
            Remove documentation references to cmd.cmds. This bypasses:
            - check to see if key is official keyword
            - converting str to python type

        """

        if isinstance(new_input, UserCommands):
            tmp_cmds = new_input.cmds
        elif isinstance(new_input, dict):
            tmp_cmds = new_input
        elif isinstance(new_input, str):
            tmp_cmds = cutils.read_config(new_input)
        else:
            raise ValueError("Cannot update with type: " + type(new_input))

        for key in tmp_cmds:
            self[key] = tmp_cmds[key]

    def validate(self):
        """
        Checks to see if all files exist and can be found by SimCADO

        Also forces any other overriding behaviour like:
            - OBS_ZENITH_DIST overrides ATMO_AIRMASS

        Returns
        -------
        valid : Bool
            True if all files are found.

        """

        # If OBS_ZENITH_DIST is specified it overrides ATMO_AIRMASS,
        if self.cmds["OBS_ZENITH_DIST"] is not None:
            airmass = zendist2airmass(self.cmds["OBS_ZENITH_DIST"])
            self.cmds["ATMO_AIRMASS"] = airmass

        # Update filenames to absolute paths
        try:
            missing_files = self._find_files()
        except ValueError:
            logging.warning("Local package database couldn't be found")
            missing_files = -1
        files_exist = len(missing_files) == 0

        # .. todo:: implement a method to check that the file formats are happy
        data_format_ok = True

        return files_exist and data_format_ok

    def writeto(self, filename="tests_commands.config"):
        """
        Write all the key-value tests_commands to an ASCII file on disk

        Parameters
        ----------
        filename : str
            file path for where the file should be saved
        """
        outstr = ""
        for group in (self.obs, self.sim,
                      self.atmo, self.scope, self.inst,
                      self.fpa, self.hxrg):
            for key in group:
                val = self[key]
                if key == "FPA_CHIP_LAYOUT" and "\n" in val:
                    val = "small"
                outstr += key.ljust(25)+"  "+str(val) + "\n"
            outstr += "\n"
        with open(filename, "w") as fd1:
            fd1.write(outstr)

    def set_instrument(self, instrument_name):
        """
        Takes the name of an INSTRUMENT package and reads the config file

        Parameters
        ----------
        instrument_name : str
            The name of a locally saved instrument package.
            Use :func:`scopesim.server.get_local_packages` to list local packages

        Examples
        --------
        ::

            >>> import scopesim as sim
            >>> cmd = sim.UserCommands()
            >>> cmd.set_instrument("MICADO")

        See Also
        --------
        :func:`scopesim.server.get_local_packages`
        :func:`scopesim.server.get_server_packages`
        :func:`scopesim.server.download_package`
        :func:`scopesim.server.list_instruments`
        :func:`scopesim.server.list_psfs`

        """

        pkg_dir = self._find_package_path(instrument_name)
        if pkg_dir is None:
            raise ValueError("{} package was not found in the folder: "
                             "\n {}".format(instrument_name, pkg_dir))

        pkg_base_file = os.path.join(pkg_dir, instrument_name+".config")
        self.update(pkg_base_file)
        self._find_files()

    def set_mode(self, mode_name):
        """
        Takes the name of an instrument MODE and reads appropriate config file

        Parameters
        ----------
        mode_name : str
            The name of a mode config file inside the instrument package
            Use ``<UserCommands>.list("modes")`` to list the available modes

        Examples
        --------
        ::

            >>> import scopesim as sim
            >>> cmd_micado = sim.UserCommands(instrument="MICADO")
            >>> cmd_micado.set_mode("MODE_MCAO_WIDE")
            >>> cmd_micado.list("modes")

        """

        if ".config" not in mode_name:
            mode_name += ".config"
        mode_path = find_file(mode_name, [self.inst_pkg_path])

        if mode_path is None:
            raise ValueError("Instrument package isn't set. mode_path is None")
        elif not os.path.exists(mode_path):
            raise ValueError("{} mode file was not found in the folder: "
                             "\n {}".format(mode_name, mode_path))

        self.update(mode_path)
        self._find_files()

    def select_filter(self, new_filter):
        """
        Sets the filter to be used by the optical train.

        Filters can be provided in several formats:
        - TransmissionCurve object
        - full file path : using the pattern "./TC_filter_<name>.dat"
        - filter name string : providing only the <name> part of the file name
            Note: this only works if an instrument config file has been read

        Acceptable filter names for an instrument package can be found by using
        ``<UserCommands>.list("filters")``

        Filter file names follow the TC_filter_<name>.dat convention and should
        consist of the two columns ``wavelength`` and ``transmission``

        Parameters
        ----------
        new_filter : str, TransmissionCurve
            The new filter : filename, acceptable filter name, or preloaded
            TransmissionCurve object

        See Also
        ``.list("filters")``

        """

        if isinstance(new_filter, sc.TransmissionCurve):
            self.cmds["INST_FILTER_TC"] = new_filter

        elif isinstance(new_filter, str):
            if os.path.exists(new_filter):
                self.cmds["INST_FILTER_TC"] = new_filter

            elif cutils.extract_filter_name(new_filter) in self.list("filters"):
                filter_name = cutils.extract_filter_name(new_filter)
                full_name = "TC_filter_{}.dat".format(filter_name)
                self.cmds["INST_FILTER_TC"] = full_name

            else:
                logging.warning("{} was not found".format(new_filter))

        else:
            raise TypeError("{} must have type `str` or "
                            "`TransmissionCurve`".format(new_filter))
        self._find_files()

    def list(self, what="modes", pattern=None):
        """
        List available choices for instrument 'modes' and 'filters'

        Parameters
        ----------
        what : str
            What aspect to list. Currently accepts:
            - "modes",
            - "filters"

        pattern : str, optional
            A specific glob pattern to search for in the instrument package
            directory: e.g. ``pattern=".TER_*.dat" `` will return the list
            of non-filter transmission curves included in the package directory

        Examples
        --------
        ::

            >>> import scopesim as sim
            >>> cmd_micado = sim.UserCommands(instrument="MICADO")
            >>> modes_list = cmd_micado.list("modes")
            >>> print(cmd_micado.list("filters"))

        """

        if "mode" in what and pattern is None:
            pattern = ".config"
        elif "filt" in what and pattern is None:
            pattern = "TC_filter_"

        if self.inst_pkg_path is None:
            logging.warning("SIM_INSTRUMENT_PACKAGE is not set")
            my_list = []
        else:
            glob_list = glob.glob(self.inst_pkg_path + "/*" + pattern + "*")
            my_list = [os.path.basename(mode).split(".")[0].replace(pattern, "")
                       for mode in glob_list]

        return my_list

    def _convert_python_types(self):
        """Convert string "none" values into python ``None`` values"""
        self.cmds = cutils.convert_dict_strings_to_python_types(self.cmds)

    def _find_files(self, silent=True):
        """
        Checks for files in the package search path

        ``_find_files`` searches the following directories in this order::

            . <i.e. local working dir>
            SIM_DATA_DIR
            SIM_INSTRUMENT_PACKAGE
            SIM_TELESCOPE_PACKAGE
            FILE_PSF_LOCAL_PATH
            SIM_SOURCE_PACKAGE

        Parameters
        ----------
        silent : bool
            If True, displays any broken paths

        Returns
        -------
        broken_paths : list

        """

        pkg_paths = ["./", self.sim_data_dir,
                     self.inst_pkg_path, self.scope_pkg_path,
                     self.psf_pkg_path, self.source_pkg_path]
        new_search_path = [pp for pp in pkg_paths if pp is not None]
        new_search_path += rc.__search_path__

        if self["SIM_INSTRUMENT_PACKAGE"] is None:
            logging.warning("Instrument package not set.")

        self["SCOPE_PSF_FILE"] = self._get_psf_path()

        broken_paths = []
        for key in self.cmds:
            if key == "OBS_OUTPUT_DIR":       # need not exist
                continue

            keyval = self.cmds[key]

            # not a string: not a filename
            if not isinstance(keyval, str):
                continue

            # ignore RC keywords that start with FILE_
            if key[:4] == "FILE":
                continue

            # If string has no extension, assume it's not a file name.
            # This is a strong assumption, but we need to guard from
            # looking for "yes", "no", "none", "scao", etc.
            # TODO Can we have a list of reserved keywords?
            if "." in keyval and len(keyval.split(".")[-1]) > 1:
                # look for the file
                fname = find_file(keyval, path=new_search_path, silent=silent)
                if fname is None:
                    broken_paths += [keyval]
                    logging.warning("Keyword "+key+" path doesn't exist: "
                                  + keyval)
                else:
                    self.cmds[key] = fname

        return broken_paths

    def _get_psf_path(self):
        """
        Find the path to the specified PSF

        Checks first if the path exists. If not, it looks in the local psf
        database for anything matching the current name in SCOPE_PSF_FILE

        """

        if self["SCOPE_PSF_FILE"] is None:
            psf_path = None
        elif ".fits" in self["SCOPE_PSF_FILE"]:
            if os.path.exists(self["SCOPE_PSF_FILE"]):
                psf_path = self["SCOPE_PSF_FILE"]
            else:
                psf_path = None
                logging.warning("PSF file not found: {}"
                              "".format(self["SCOPE_PSF_FILE"]))
        elif svr.get_path(self["SCOPE_PSF_FILE"]) is not None:
            psf_path = svr.get_path(self["SCOPE_PSF_FILE"])
        else:
            logging.warning("PSF file not found in the DB: {}"
                          "".format(self["SCOPE_PSF_FILE"]))
            psf_path = None

        return psf_path

    def _update_lam_extremes(self):
        """Gets the current wave_min and wave_max"""
        # if SIM_USE_FILTER_LAM is true, then use the filter curve to set the
        # wavelength boundaries where the filter is < SIM_FILTER_THRESHOLD

        if self.cmds["SIM_USE_FILTER_LAM"].lower() == "yes":
            if isinstance(self.cmds["INST_FILTER_TC"], str):
                tc_filt = sc.TransmissionCurve(
                    filename=find_file(self.cmds['INST_FILTER_TC']))
            else:
                tc_filt = self.cmds["INST_FILTER_TC"]
            mask = np.where(tc_filt.val > self.cmds["SIM_FILTER_THRESHOLD"])[0]
            imin = np.max((mask[0] - 1, 0))
            imax = np.min((mask[-1] + 1, len(tc_filt.lam) - 1))
            lam_min, lam_max = tc_filt.lam[imin], tc_filt.lam[imax]
        else:
            lam_min = self.cmds["SIM_LAM_MIN"]
            lam_max = self.cmds["SIM_LAM_MAX"]

        return lam_min, lam_max

    def _get_lam_bin_edges(self, lam_min, lam_max):
        """
        Generates an array with the bin edges of the layers in spectral space

        Parameters
        ----------
        lam_min, lam_max : float
            [um] the minimum and maximum wavelengths of the filter range

        Notes
        -------
        Atmospheric diffraction causes blurring in an image. To model this
        effect the spectra from a ``Source`` object are cut into bins based on
        how far the photons at different wavelength are diffracted from the
        image center. The threshold for defining a new layer based on the how
        far a certain bin will move is given by ``SIM_ADC_SHIFT_THRESHOLD``. The
        default value is 1 pixel.

        The PSF also causes blurring as it spreads out over a bandpass. This
        also needed to be taken into account
        """

        if self.cmds["SIM_VERBOSE"] == "yes":
            print("Determining lam_bin_edges")

        effectiveness = self.cmds["INST_ADC_PERFORMANCE"] / 100.

        # This is redundant because also need to look at the PSF width
        # if effectiveness == 1.:
        #    lam_bin_edges = np.array([lam_min, lam_max])
        #    return lam_bin_edges

        shift_threshold = self.cmds["SIM_ADC_SHIFT_THRESHOLD"]

        # get the angle shift for each slice
        lam = np.arange(lam_min, lam_max + 1E-7, 0.001)
        zenith_distance = airmass2zendist(self.cmds["ATMO_AIRMASS"])
        angle_shift = atmospheric_refraction(lam,
                                             zenith_distance,
                                             self.cmds["ATMO_TEMPERATURE"],
                                             self.cmds["ATMO_REL_HUMIDITY"],
                                             self.cmds["ATMO_PRESSURE"],
                                             self.cmds["SCOPE_LATITUDE"],
                                             self.cmds["SCOPE_ALTITUDE"])

        # convert angle shift into number of pixels
        # pixel shifts are defined with respect to last slice
        rel_shift = (angle_shift - angle_shift[-1]) / self.pix_res
        rel_shift *= (1. - effectiveness)
        if np.max(np.abs(rel_shift)) > 1000:
            raise ValueError("Pixel shifts too great (>1000), check units")

        # Rotate by the paralytic angle
        int_shift = np.array(rel_shift / shift_threshold, dtype=np.int)
        idx = [np.where(int_shift == i)[0][0]
               for i in np.unique(int_shift)[::-1]]
        lam_bin_edges_adc = np.array(lam[idx + [len(lam)-1]])

        # Now check to see if the PSF blurring is the controlling factor. If so,
        # take the lam_bin_edges for the PSF blurring

        diam = self.diameter
        d_ang = self.pix_res * shift_threshold

        # .. todo:: get rid of hard coded diameter of MICADO FOV
        diam_arcsec = 1.22 * 53 * 3600
        lam_bin_edges_psf = [lam_min]
        ang0 = (lam_min*1E-6) / diam * diam_arcsec

        i = 1
        while lam_bin_edges_psf[-1] < lam_max:
            lam_bin_edges_psf += [(ang0 + d_ang*i) * diam / diam_arcsec * 1E6]
            i += 1
            if i > 1000:
                raise ValueError("lam_bin_edges needs >1000 values")
        lam_bin_edges_psf[-1] = lam_max

        lam_bin_edges = np.unique(np.concatenate(
            (np.round(lam_bin_edges_psf, 3),
             np.round(lam_bin_edges_adc, 3))))

        if self.cmds["SIM_VERBOSE"] == "yes":
            print("PSF edges were", np.round(lam_bin_edges_psf, 3))
            print("ADC edges were", np.round(lam_bin_edges_adc, 3))
            print("All edges were", np.round(lam_bin_edges, 3))

        return lam_bin_edges

    def _find_package_path(self, pkg_name):
        """Looks in local package DB for path to the package data"""

        pkg_entry = svr.find_package_on_disk(pkg_name)
        if pkg_entry is None:
            full_pkg_path = None
        else:
            pkg_path = pkg_entry["path"].replace(".zip", "")
            downloads_path = self["FILE_LOCAL_DOWNLOADS_PATH"]
            full_pkg_path = os.path.join(downloads_path, pkg_path)

            if not os.path.exists(full_pkg_path):
                raise ValueError("{} is given in the database, but doesn't "
                                 "exist".format(full_pkg_path))

        return full_pkg_path

    def __str__(self):
        if self.cmds["CONFIG_USER"] is not None:
            return "A dictionary of tests_commands compiled from " + \
                                                        self.cmds["CONFIG_USER"]
        return "A dictionary of default tests_commands"

    def __iter__(self):
        return self.cmds.__iter__()

    def __next__(self):
        return self.cmds.__next__()

    def __getitem__(self, key):
        if cutils.is_item_subcategory(key, self.cmds):
            return cutils.get_subcategory(key, self.cmds)
        else:
            return self.cmds[key]

    def __setitem__(self, key, val):
        if key not in self.cmds:
            logging.warning("{} not in self.keys. Ignoring.".format(key))
            return None

        self.cmds[key] = cutils.str_to_python_type(val)

    # # Doesn't work because of a recursion error with copy.deepcopy
    # def __getattr__(self, item):
    #     """For things like subcategories, e.g. self.atmo """
    #     if item in self.cmds:
    #         attr = self.cmds[item]
    #     elif cutils.is_item_subcategory(item, self.cmds):
    #         attr = cutils.get_subcategory(item, self.cmds)
    #     else:
    #         logging.warning("{} doesn't exist. Capital letters?". format(item))
    #         raise AttributeError(item)
    #     return attr

    @property
    def yaml_docs(self):
        sys_descrip = self.cmds["SIM_SYSTEM_DESCRIPTION"]
        with open(find_file(sys_descrip)) as filename:
            yaml_files = yaml.full_load(filename)
        for element in yaml_files:
            with open(find_file(yaml_files[element])) as filename:
                yaml_files[element] = yaml.full_load(filename)

        return yaml_files

    @property
    def keys(self):
        """Return the keys in the ``UserCommands.cmds`` dictionary"""
        return self.cmds.keys()

    @property
    def values(self):
        """Return the values in the ``UserCommands.cmds`` dictionary"""
        return self.cmds.values()

    # All the convenience attributes from the update_attributes function
    @property
    def diameter(self):
        if self.mirrors_telescope is not None:
            i = np.where(self.mirrors_telescope["Mirror"] == "M1")[0][0]
            return self.mirrors_telescope["Outer"][i]

    @property
    def area(self):
        if self.mirrors_telescope is not None:
            i = np.where(self.mirrors_telescope["Mirror"] == "M1")[0][0]
            return np.pi / 4 * (self.mirrors_telescope["Inner"][i] ** 2 -
                                self.mirrors_telescope["Inner"][i] ** 2)

    @property
    def total_wfe(self):
        if self.cmds["INST_WFE"] is not None:
            if isinstance(self.cmds["INST_WFE"], str):
                wfe_list = ioascii.read(self.cmds["INST_WFE"])
                wfe = wfe_list["wfe_rms"]
                num = wfe_list["no_surfaces"]
            elif isinstance(self.cmds["INST_WFE"], (float, int)):
                wfe, num = float(self.cmds["INST_WFE"]), 1
            print(num, wfe)
            tot_wfe = np.sqrt(np.sum(num * wfe**2))
        else:
            logging.warning("INST_WFE is None. Returning zero wavefront error")
            tot_wfe = 0

        return tot_wfe

    @property
    def mirrors_telescope(self):
        if self.cmds["SCOPE_MIRROR_LIST"] is not None:
            return ioascii.read(self.cmds["SCOPE_MIRROR_LIST"])

    @property
    def mirrors_ao(self):
        if self.cmds["INST_MIRROR_AO_LIST"] is not None:
            return ioascii.read(self.cmds["INST_MIRROR_AO_LIST"])

    @property
    def fpa_res(self):
        return self.cmds["SIM_PIXEL_SCALE"]

    @property
    def pix_res(self):
        return self.fpa_res / self.cmds["SIM_OVERSAMPLING"]

    @property
    def lam_res(self):
        return self.cmds["SIM_LAM_TC_BIN_WIDTH"]

    @property
    def lam(self):
        lam_min, lam_max = self._update_lam_extremes()
        return np.arange(lam_min, lam_max + 1E-7, self.lam_res)

    @property
    def _lam_min(self):
        return self.lam[0]

    @property
    def _lam_max(self):
        return self.lam[-1]

    @property
    def lam_bin_edges(self):
        lam_min, lam_max = self._update_lam_extremes()
        return self._get_lam_bin_edges(lam_min, lam_max)

    @property
    def lam_bin_centers(self):
        return 0.5 * (self.lam_bin_edges[1:] + self.lam_bin_edges[:-1])

    @property
    def exptime(self):
        return self.cmds["OBS_EXPTIME"]

    @property
    def verbose(self):
        return self.cmds["SIM_VERBOSE"]

    @property
    def sim_data_dir(self):
        return self.cmds["SIM_DATA_DIR"]

    @property
    def inst_pkg_path(self):
        return self._find_package_path(self["SIM_INSTRUMENT_PACKAGE"])

    @property
    def scope_pkg_path(self):
        return self._find_package_path(self["SIM_TELESCOPE_PACKAGE"])

    @property
    def psf_pkg_path(self):
        return self._find_package_path(self["FILE_PSF_LOCAL_PATH"])

    @property
    def source_pkg_path(self):
        return self._find_package_path(self["SIM_SOURCE_PACKAGE"])

    def _default_data(self):
        # .. todo:: The only thing that was missing was custom detector layouts
        # So implement a custom detector read out function in detector, not here
        pass
