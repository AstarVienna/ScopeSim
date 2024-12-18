# -*- coding: utf-8 -*-

import copy

from datetime import datetime

import numpy as np
from scipy.interpolate import interp1d
from astropy import units as u
from astropy.wcs import WCS

from tqdm.auto import tqdm

from synphot import SourceSpectrum, Empirical1D
from synphot.units import PHOTLAM

from astar_utils.nested_mapping import RecursiveNestedMapping, is_bangkey, recursive_update

from .optics_manager import OpticsManager
from .fov_manager import FOVManager
from .image_plane import ImagePlane
from ..commands.user_commands import UserCommands
from ..detector import DetectorManager
from ..effects import ExtraFitsKeywords
from ..utils import from_currsys, top_level_catch, get_logger
from ..source.source_templates import empty_sky
from .. import __version__


logger = get_logger(__name__)

import multiprocessing as mp

N_PROCESSES = mp.cpu_count() - 1
USE_MULTIPROCESSING = False

def extract_source(fov, source):
    fov.extract_from(source)
    return fov

def view_fov(fov, hdu_type):
    fov.view(hdu_type)
    return fov

def apply_fov_effects(fov, fov_effects):
    for effect in fov_effects:
        fov = effect.apply_to(fov)
    return fov


class OpticalTrain:
    """
    The main class for controlling a simulation.

    Parameters
    ----------
    cmds : UserCommands, str
        If the name of an instrument is passed, OpticalTrain tries to find the
        instrument package, and internally creates the UserCommands object

    Examples
    --------
    Create an optical train::

        >>> import scopesim as im
        >>> cmd = sim.UserCommands("MICADO")
        >>> opt = sim.OpticalTrain(cmd)

    Observe a Source object::

        >>> src = sim.source.source_templates.empty_sky()
        >>> opt.observe(src)
        >>> hdus = opt.readout()

    If no Source is specified, an empty field is observe, so that the
    following is equivalent to the commands above::
        >>> opt.observe()
        >>> hdus = opt.readout()

    List the effects modelled in an OpticalTrain::

        >>> print(opt.effects)

    Effects can be accessed by using the name of the effect::

        >>> print(opt["dark_current"])

    To include or exclude an effect during a simulation run, use the
    ``.include`` attribute of the effect::

         >>> opt["dark_current"].include = False

    Data used by an Effect object is contained in the ``.data`` attribute,
    while other information is contained in the ``.meta`` attribute::

        >>> opt["dark_current"].data
        >>> opt["dark_current"].meta

    Meta data values can be set by either using the ``.meta`` attribute
    directly::

        >>> opt["dark_current"].meta["value"] = 0.5

    or by passing a dictionary (with one or multiple entries) to the
    OpticalTrain object::

        >>> opt["dark_current"] = {"value": 0.75, "dit": 30}

    """

    @top_level_catch
    def __init__(self, cmds=None):
        self.cmds = cmds
        self._description = self.__repr__()
        self.optics_manager = None
        self.fov_manager = None
        self.image_planes = []
        self.detector_managers = []
        self.yaml_dicts = None
        self._last_source = None

        if cmds is not None:
            self.load(cmds)

    def load(self, user_commands):
        """
        (Re)Load an OpticalTrain with a new set of UserCommands.

        Parameters
        ----------
        user_commands : UserCommands or str

        """
        if isinstance(user_commands, str):
            user_commands = UserCommands(use_instrument=user_commands)
        elif isinstance(user_commands, UserCommands):
            user_commands = copy.deepcopy(user_commands)
        else:
            raise ValueError("user_commands must be a UserCommands or str object "
                             f"but is {type(user_commands)}")

        # HACK: The clean solution of new_child doesn't work because all the
        #       extras of UserCommands like yaml_dicts and so on are not
        #       passed to the new instance. Ideally fix that at some point...
        # self.cmds = user_commands.new_child(RecursiveNestedMapping(title="CurrObs"))
        self.cmds = user_commands
        self.cmds.maps = [RecursiveNestedMapping(title="CurrObs"), *self.cmds.maps]
        # FIXME: Setting rc.__currsys__ to user_commands causes many problems:
        #        UserCommands used NestedMapping internally, but is itself not
        #        an instance or subclas thereof. So rc.__currsys__ actually
        #        changes type as a result of this line. On one hand, some other
        #        code relies on this change, i.e. uses attributes from
        #        UserCommands via rc.__currsys__, but on the other hand some
        #        tests (now with proper patching) fail because of this type
        #        change. THIS IS A PROBLEM!
        # NOTE: All tests pass without setting rc.__currsys__ to user_commands.
        #       Nevertheless, I'm a bit reluctant to removing this code just
        #       yet. So it is commented out.
        # rc.__currsys__ = user_commands
        self.yaml_dicts = self.cmds.yaml_dicts
        self.optics_manager = OpticsManager(self.yaml_dicts, self.cmds)
        self.update()

    def update(self, **kwargs):
        """
        Update the user-defined parameters and remake main internal classes.

        Parameters
        ----------
        kwargs : expanded dict
            Any keyword-value pairs from a config file

        """
        self.optics_manager.update(**kwargs)
        opt_man = self.optics_manager

        self.fov_manager = FOVManager(opt_man.fov_setup_effects, cmds=self.cmds,
                                      **kwargs)
        self.image_planes = [ImagePlane(hdr, self.cmds, **kwargs)
                             for hdr in opt_man.image_plane_headers]
        self.detector_managers = [DetectorManager(det_list, cmds=self.cmds, **kwargs)
                                for det_list in opt_man.detector_setup_effects]

        # Move everything from CurrObs to CurrSys, so CurrObs is clean for
        # .observe and .readout. This is necessary because the setup and
        # potentially update need to be able to modify CurrSys, but CurrObs
        # needs to be already in place when creating e.g. effects, but then
        # writing to cmds will always use CurrObs. The alternative is manually
        # writing to cmds.maps[1] everywhere during the setup, but these two
        # lines make sure everything does indeed go into CurrSys.
        # self.cmds.maps[1].update(self.cmds.maps[0])
        # self.cmds.maps[0].clear()
        # HACK: recursive_update is needed to avoid overwriting emtpy !ABCs
        # TODO: or is it??
        self.cmds.maps[1].dic = recursive_update(self.cmds.maps[1].dic, self.cmds.maps[0].dic)
        self.cmds.maps[0].dic.clear()

    @top_level_catch
    def observe(self, orig_source=None, update=True, **kwargs):
        """
        Main controlling method for observing ``Source`` objects.

        Parameters
        ----------
        orig_source : Source
        update : bool
            Reload optical system
        kwargs : expanded dict
            Any keyword-value pairs from a config file

        Notes
        -----
        When orig_source is None (e.g. when writing opt.observe()), an
        empty field is observed (internally created with empty_sky()).

        How the list of Effects is split between the 5 main tasks:

        - Make a FOV list - z_order = 0..99
        - Make a image plane - z_order = 100..199
        - Apply Source altering effects - z_order = 200..299
        - Apply FOV specific (3D) effects - z_order = 300..399
        - Apply FOV-independent (2D) effects - z_order = 400..499
        - [Apply detector plane (0D, 2D) effects - z_order = 500..599]

        .. todo:: List is out of date - update

        """
        if update:
            self.update(**kwargs)

        # self.set_focus(**kwargs)    # put focus back on current instrument package

        # Make a copy of the Source and prepare for observation (convert to
        # internally used units, sample to internal wavelength grid)
        if orig_source is None:
            source = empty_sky()
            logger.info("Observing empty field")
        else:
            source = orig_source.make_copy()
        source = self.prepare_source(source)

        # [1D - transmission curves]
        for effect in self.optics_manager.source_effects:
            source = effect.apply_to(source)

        # [3D - Atmospheric shifts, PSF, NCPAs, Grating shift/distortion]

        # START OF MULTIPROCESSING
        if USE_MULTIPROCESSING:

            fovs = self.fov_manager.fovs
            fov_effects = self.optics_manager.fov_effects
            hdu_type = "cube" if self.fov_manager.is_spectroscope else "image"

            with mp.Pool(processes=N_PROCESSES) as pool:
                fovs = pool.starmap(extract_source,
                                    zip(fovs, [source] * len(fovs)))

            with mp.Pool(processes=N_PROCESSES) as pool:
                fovs = pool.starmap(view_fov,
                                    zip(fovs, [hdu_type] * len(fovs)))

            with mp.Pool(processes=N_PROCESSES) as pool:
                fovs = pool.starmap(apply_fov_effects,
                                    zip(fovs, [fov_effects] * len(fovs)))

        # OLD SINGLE CORE CODE
        else:

            fovs = self.fov_manager.fovs
            nobar = len(fovs) <= 1
            for fov in tqdm(fovs, desc=" FOVs", position=0, disable=nobar):
                # print("FOV", fov_i+1, "of", n_fovs, flush=True)
                # .. todo: possible bug with bg flux not using plate_scale
                #          see fov_utils.combine_imagehdu_fields
                fov.extract_from(source)

                hdu_type = "cube" if self.fov_manager.is_spectroscope else "image"
                fov.view(hdu_type)
                foveffs = self.optics_manager.fov_effects
                nobar = len(foveffs) <= 1
                for effect in tqdm(foveffs, disable=nobar,
                                   desc=" FOV effects", position=1):#, leave=False):
                    fov = effect.apply_to(fov)

                fov.flatten()
                self.image_planes[fov.image_plane_id].add(fov.hdu, wcs_suffix="D")
                # ..todo: finish off the multiple image plane stuff

        # END OF MULTIPROCESSING

        # [2D - Vibration, flat fielding, chopping+nodding]
        impeffs = self.optics_manager.image_plane_effects
        nobar = len(impeffs) <= 1
        for effect in tqdm(impeffs, disable=nobar,
                           desc=" Image Plane effects"):
            for ii, image_plane in enumerate(self.image_planes):
                self.image_planes[ii] = effect.apply_to(image_plane)

        self._last_fovs = fovs
        self._last_source = source

    def prepare_source(self, source):
        """
        Prepare source for observation.

        The method is currently applied to cube fields only.
        The source data are converted to internally used units (PHOTLAM).
        The source data are interpolated to the waveset used by the FieldOfView
        This is necessary when the source data are sampled on a coarser grid
        than used internally, or if the source data are sampled on irregular
        wavelengths.
        For cube fields, the method assumes that the wavelengths at which the
        cube is sampled is provided explicitely as attribute `wave` if the cube
        ImageHDU.
        """
        # Convert to PHOTLAM per arcsec2
        # TODO: this is not sufficiently general
        # TODO: Maybe move this to source_fields??

        for ispec, spec in source.spectra.items():
            # Put on fov wavegrid
            wave_min = min(fov.meta["wave_min"] for fov in self.fov_manager.fovs)
            wave_max = max(fov.meta["wave_max"] for fov in self.fov_manager.fovs)
            wave_unit = u.Unit(from_currsys("!SIM.spectral.wave_unit", self.cmds))
            dwave = from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)  # Not a quantity
            fov_waveset = np.arange(wave_min.value, wave_max.value, dwave) * wave_unit
            fov_waveset = fov_waveset.to(u.um)

            source.spectra[ispec] = SourceSpectrum(Empirical1D,
                                                   points=fov_waveset,
                                                   lookup_table=spec(fov_waveset))

        for cube in source.cube_fields:
            header, data, wave = cube.header, cube.data, cube.wave

            # Need to check whether BUNIT is per arcsec2 or per pixel
            inunit = u.Unit(header["BUNIT"])
            data = data.astype(np.float32) * inunit
            factor = 1
            for base, power in zip(inunit.bases, inunit.powers):
                if (base**power).is_equivalent(u.arcsec**(-2)):
                    conversion = (base**power).to(u.arcsec**(-2)) / base**power
                    data *= conversion
                    factor = u.arcsec**(-2)

            data = data.to(PHOTLAM,
                           equivalencies=u.spectral_density(wave[:, None, None]))

            if factor == 1:    # Normalise to 1 arcsec2 if not a spatial density
                # ..todo: lower needed because "DEG" is not understood, this is ugly
                pixarea = (header["CDELT1"] * u.Unit(header["CUNIT1"].lower()) *
                           header["CDELT2"] * u.Unit(header["CUNIT2"].lower())).to(u.arcsec**2)
                data = data / pixarea.value    # cube is per arcsec2

            data = (data * factor).value

            cube.header["BUNIT"] = "PHOTLAM/arcsec2"    # ..todo: make this more explicit?

            # The imageplane_utils like to have the spatial WCS in units of "deg". Ensure
            # that the cube is passed on accordingly
            cube.header["CDELT1"] = header["CDELT1"] * u.Unit(header["CUNIT1"].lower()).to(u.deg)
            cube.header["CDELT2"] = header["CDELT2"] * u.Unit(header["CUNIT2"].lower()).to(u.deg)
            cube.header["CUNIT1"] = "deg"
            cube.header["CUNIT2"] = "deg"

            # Put on fov wavegrid
            wave_min = min(fov.meta["wave_min"] for fov in self.fov_manager.fovs)
            wave_max = max(fov.meta["wave_max"] for fov in self.fov_manager.fovs)
            wave_unit = u.Unit(from_currsys("!SIM.spectral.wave_unit", self.cmds))
            dwave = from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)  # Not a quantity
            fov_waveset = np.arange(wave_min.value, wave_max.value, dwave) * wave_unit
            fov_waveset = fov_waveset.to(u.um)

            # Interpolate into new data cube.
            # This is done layer by layer for memory reasons.
            new_data = np.zeros((fov_waveset.shape[0], data.shape[1], data.shape[2]),
                                dtype=np.float32)
            for j in range(data.shape[1]):
                cube_interp = interp1d(wave.to(u.um).value, data[:, j, :],
                                       axis=0, kind="linear",
                                       bounds_error=False, fill_value=0)
                new_data[:, j, :] = cube_interp(fov_waveset.value)

            cube.data = new_data
            cube.header["CTYPE3"] = "WAVE"
            cube.header["CRPIX3"] = 1
            cube.header["CRVAL3"] = wave_min.value
            cube.header["CDELT3"] = dwave
            cube.header["CUNIT3"] = wave_unit.name

            cube.wcs = WCS(cube.field)

        return source

    @top_level_catch
    def readout(self, filename=None, reset=True, **kwargs):
        """
        Produce detector readouts for the observed image.

        Parameters
        ----------
        filename : str, optional
            Where to save the FITS file
        kwargs

        Returns
        -------
        hdu : fits.HDUList

        Notes
        -----
        - Apply detector plane (0D, 2D) effects - z_order = 500..599

        """
        if reset:
            self.cmds.clear()

        # Hack to make sure AutoExposure and Quantization work properly.
        # Should probably be removed once #428 is fixed properly.
        # if kwargs.get('exptime', None) is not None:
        #     # self.cmds.pop("!OBS.dit", None)
        #     # self.cmds.pop("!OBS.ndit", None)
        #     self.cmds["!OBS.dit"] = None
        #     self.cmds["!OBS.ndit"] = None
        # TODO: This is still hacky but seems to work for now...
        params = {"exptime": None}
        params.update(kwargs)
        if params["exptime"] is not None and params.get("dit") is None:
            params.update(
                dit=None,
                ndit=None,
            )
        else:
            params.update(
                exptime=params.get("exptime") or self.cmds.get("!OBS.exptime"),
                dit=params.get("dit") or self.cmds.get("!OBS.dit"),
                ndit=params.get("ndit") or self.cmds.get("!OBS.ndit"),
            )

        for key, value in params.items():
            if is_bangkey(key):
                # skip if already explicitly !ABC.xyz
                logger.debug("%s: %s given in kwargs", key, value)
                continue
            # otherwise fallback to !OBS
            logger.debug("%s: %s given in kwargs, put in !OBS", key, value)
            self.cmds[f"!OBS.{key}"] = value

        hduls = []
        for i, detector_array in enumerate(self.detector_managers):
            array_effects = self.optics_manager.detector_array_effects
            dtcr_effects = self.optics_manager.detector_effects
            hdul = detector_array.readout(
                self.image_planes, array_effects, dtcr_effects)

            fits_effects = self.optics_manager.get_all(ExtraFitsKeywords)
            if len(fits_effects) > 0:
                for effect in fits_effects:
                    hdul = effect.apply_to(hdul, optical_train=self)
            else:
                try:
                    hdul = self.write_header(hdul)
                except Exception:
                    logger.exception("Header update failed, data will be "
                                     "saved with incomplete header. See stack "
                                     "trace for details.")

            if filename is not None and isinstance(filename, str):
                fname = filename
                if len(self.detector_managers) > 1:
                    fname = f"{i}_{filename}"
                hdul.writeto(fname, overwrite=True)

            hduls.append(hdul)

        # Create a copy of the readouts, because they are still referenced here.
        # That is, self.detector_managers[0][0].hdu is hdul[1].
        # Subsequent readouts would therefor overwrite earlier readouts
        # without the copy.
        # It is necessary to keep the detector readouts in memory, because
        # this makes it possible to add support for up-the-ramp sampling,
        # that is, subsequent readouts without reseting the detector inbetween.
        return copy.deepcopy(hduls)

    def write_header(self, hdulist):
        """Write meaningful header to simulation product."""
        # Primary hdu
        pheader = hdulist[0].header
        pheader["DATE"] = datetime.now().isoformat(timespec="seconds")
        pheader["ORIGIN"] = f"Scopesim {__version__}"
        pheader["INSTRUME"] = from_currsys("!OBS.instrument", self.cmds)
        pheader["INSTMODE"] = ", ".join(from_currsys("!OBS.modes", self.cmds))
        pheader["TELESCOP"] = from_currsys("!TEL.telescope", self.cmds)
        pheader["LOCATION"] = from_currsys("!ATMO.location", self.cmds)

        # Source information taken from first only.
        # ..todo: What if source is a composite?
        srcfield = self._last_source.fields[0]
        if type(srcfield).__name__ == "Table":
            pheader["SOURCE"] = "Table"
        elif type(srcfield).__name__ == "ImageHDU":
            if "BG_SURF" in srcfield.header:
                pheader["SOURCE"] = srcfield.header["BG_SURF"]
            else:
                try:
                    pheader["SOURCE"] = srcfield.header["FILENAME"]
                except KeyError:
                    pheader["SOURCE"] = "ImageHDU"

        # Image hdul
        # ..todo: currently only one, update for detector arrays
        # ..todo: normalise filenames - some need from_currsys, some need Path(...).name
        #         this should go into a function so as to reduce clutter here.
        iheader = hdulist[1].header
        iheader["EXPTIME"] = from_currsys("!OBS.exptime", self.cmds), "[s]"
        iheader["DIT"] = from_currsys("!OBS.dit", self.cmds), "[s]"
        iheader["NDIT"] = from_currsys("!OBS.ndit", self.cmds)
        iheader["BUNIT"] = "e", "per EXPTIME"
        iheader["PIXSCALE"] = from_currsys("!INST.pixel_scale", self.cmds), "[arcsec]"

        for eff in self.optics_manager.detector_setup_effects:
            efftype = type(eff).__name__

            if efftype == "DetectorList" and eff.include:
                iheader["DETECTOR"] = eff.meta["detector"]

        for eff in self.optics_manager.detector_array_effects:
            efftype = type(eff).__name__

            if (efftype == "DetectorModePropertiesSetter" and
                eff.include):
                # ..todo: can we write this into currsys?
                iheader["DET_MODE"] = (eff.meta["detector_readout_mode"],
                                       "detector readout mode")
                iheader["MINDIT"] = from_currsys("!DET.mindit", self.cmds), "[s]"
                iheader["FULLWELL"] = from_currsys("!DET.full_well", self.cmds), "[s]"
                iheader["RON"] = from_currsys("!DET.readout_noise", self.cmds), "[e]"
                iheader["DARK"] = from_currsys("!DET.dark_current", self.cmds), "[e/s]"

        ifilter = 1   # Counts filter wheels
        isurface = 1  # Counts surface lists
        for eff in self.optics_manager.source_effects:
            efftype = type(eff).__name__

            if efftype == "ADCWheel" and eff.include:
                iheader["ADC"] = eff.current_adc.meta["name"]

            if efftype == "FilterWheel" and eff.include:
                iheader[f"FILTER{ifilter}"] = (eff.current_filter.meta["name"],
                                               eff.meta["name"])
                ifilter += 1

            if efftype == "SlitWheel" and eff.include:
                iheader["SLIT"] = (eff.current_slit.meta["name"],
                                   eff.meta["name"])

            if efftype == "PupilTransmission" and eff.include:
                iheader["PUPTRANS"] = (from_currsys("!OBS.pupil_transmission", self.cmds),
                                       "cold stop, pupil transmission")

            if efftype == "SkycalcTERCurve" and eff.include:
                iheader["ATMOSPHE"] = "Skycalc", "atmosphere model"
                iheader["LOCATION"] = eff.meta["location"]
                iheader["AIRMASS"] = eff.meta["airmass"]
                iheader["TEMPERAT"] = eff.meta["temperature"], "[degC]"
                iheader["HUMIDITY"] = eff.meta["humidity"]
                iheader["PRESSURE"] = eff.meta["pressure"], "[hPa]"
                iheader["PWV"] = eff.meta["pwv"], "precipitable water vapour"

            if efftype == "SurfaceList" and eff.include:
                iheader[f"SURFACE{isurface}"] = eff.meta["name"]
                isurface += 1

        return hdulist

    def shutdown(self):
        """
        Shut down the instrument.

        This method closes all open file handles and should be called when the
        optical train is no longer needed.
        """
        for effect_name in self.effects["name"]:
            try:
                self[effect_name]._file.close()
            except AttributeError:
                pass

        self._description = "The instrument has been shut down."


    @property
    def effects(self):
        return self.optics_manager.list_effects()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.cmds!r})"

    def __str__(self):
        return self._description

    def _repr_pretty_(self, p, cycle):
        """For ipython."""
        if cycle:
            p.text(f"{self.__class__.__name__}(...)")
        else:
            p.text(f"{self.__class__.__name__} ")
            p.text(f"for {self.cmds['!OBS.instrument']} ")
            p.text(f"@ {self.cmds['!TEL.telescope']}:")
            p.breakable()
            p.text("UserCommands:")
            p.breakable()
            p.pretty(self.cmds)
            p.breakable()
            p.text("OpticalElements:")
            with p.indent(2):
                for item in self:
                    p.breakable()
                    p.pretty(item)
            p.breakable()
            p.text("DetectorManagers:")
            with p.indent(2):
                for item in self.detector_managers:
                    p.breakable()
                    p.pretty(item)
            p.breakable()
            p.text("Effects:")
            p.breakable()
            with p.indent(2):
                p.pretty(self.effects)

    def __getitem__(self, item):
        return self.optics_manager[item]

    def __setitem__(self, key, value):
        self.optics_manager[key] = value

    def __contains__(self, key):
        return key in self.optics_manager

    def report(self):
        pass
        # user commands report
        #   package dependencies
        #   modes names
        #   default modes
        #   yaml hierarchy
        #
        # optics_manager
        #   derived properties
        #   system transmission curve
        #   list of effects
        #
        # etc
        #   limiting magnitudes
        #
