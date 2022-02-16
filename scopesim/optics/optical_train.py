from copy import deepcopy
from shutil import copyfileobj
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

from synphot.units import PHOTLAM

from .optics_manager import OpticsManager
from .fov_manager import FOVManager
from .image_plane import ImagePlane
from ..commands.user_commands import UserCommands
from ..detector import DetectorArray
from ..source.source import Source
from ..utils import from_currsys
from .. import rc
from . import fov_utils as fu


class OpticalTrain:
    """
    The main class for controlling a simulation

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

    List the effects modelled in an OpticalTrain::

        >>> print(opt.effects)

    Effects can be accessed by using the name of the effect::

        >>> print(opt["dark_current"])

    To include or exclude an effect during a simulation run, use the
    ``.include`` attribute of the effect::

         >>> opt["dark_current"].include = False

    Data used by an Effect object is contained in the ``.data`` attribute, while
    other information is contained in the ``.meta`` attribute::

        >>> opt["dark_current"].data
        >>> opt["dark_current"].meta

    Meta data values can be set by either using the ``.meta`` attribute
    directly::

        >>> opt["dark_current"].meta["value"] = 0.5

    or by passing a dictionary (with one or multiple entries) to the
    OpticalTrain object::

        >>> opt["dark_current"] = {"value": 0.75, "dit": 30}

    """
    def __init__(self, cmds=None):

        self._description = self.__repr__()
        self.cmds = cmds
        self.optics_manager = None
        self.fov_manager = None
        self.image_planes = []
        self.detector_arrays = []
        self.yaml_dicts = None
        self._last_source = None


        if cmds is not None:
            self.load(cmds)

    def load(self, user_commands):
        """
        (Re)Loads an OpticalTrain with a new set of UserCommands

        Parameters
        ----------
        user_commands : UserCommands

        """

        if isinstance(user_commands, str):
            user_commands = UserCommands(use_instrument=user_commands)

        if not isinstance(user_commands, UserCommands):
            raise ValueError("user_commands must be a UserCommands object: "
                             "{}".format(type(user_commands)))

        self.cmds = user_commands
        rc.__currsys__ = user_commands
        self.yaml_dicts = rc.__currsys__.yaml_dicts
        self.optics_manager = OpticsManager(self.yaml_dicts)
        self.update()

    def update(self, **kwargs):
        """
        Update the user-defined parameters and remake the main internal classes

        Parameters
        ----------
        kwargs : expanded dict
            Any keyword-value pairs from a config file

        """
        self.optics_manager.update(**kwargs)
        opt_man = self.optics_manager

        self.fov_manager = FOVManager(opt_man.fov_setup_effects, **kwargs)
        self.image_planes = [ImagePlane(hdr, **kwargs)
                             for hdr in opt_man.image_plane_headers]
        self.detector_arrays = [DetectorArray(det_list, **kwargs)
                                for det_list in opt_man.detector_setup_effects]


    def observe(self, orig_source, update=True, **kwargs):
        """
        Main controlling method for observing ``Source`` objects

        Parameters
        ----------
        orig_source : Source
        update : bool
            Reload optical system
        kwargs : expanded dict
            Any keyword-value pairs from a config file

        Notes
        -----
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

        self.set_focus(**kwargs)    # put focus back on current instrument package

        # Make a copy of the Source and prepare for observation (convert to
        # internally used units, sample to internal wavelength grid)
        source = orig_source.make_copy()
        source = self.prepare_source(source)

        # [1D - transmission curves]
        for effect in self.optics_manager.source_effects:
            source = effect.apply_to(source)

        # [3D - Atmospheric shifts, PSF, NCPAs, Grating shift/distortion]
        fovs = self.fov_manager.fovs
        for fov in fovs:
            # print("FOV", fov_i+1, "of", n_fovs, flush=True)
            # .. todo: possible bug with bg flux not using plate_scale
            #          see fov_utils.combine_imagehdu_fields
            fov.extract_from(source)

            hdu_type = "cube" if self.fov_manager.is_spectroscope else "image"
            fov.view(hdu_type)
            for effect in self.optics_manager.fov_effects:
                fov = effect.apply_to(fov)

            fov.flatten()
            self.image_planes[fov.image_plane_id].add(fov.hdu, wcs_suffix="D")
            # ..todo: finish off the multiple image plane stuff

        # [2D - Vibration, flat fielding, chopping+nodding]
        for effect in self.optics_manager.image_plane_effects:
            for ii in range(len(self.image_planes)):
                self.image_planes[ii] = effect.apply_to(self.image_planes[ii])

        self._last_source = source


    def prepare_source(self, source):
        """
        Prepare source for observation

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
        # ..todo: this is not sufficiently general

        for cube in source.cube_fields:
            header, data, wave = cube.header, cube.data, cube.wave

            # Need to check whether BUNIT is per arcsec2 or per pixel
            inunit = u.Unit(header['BUNIT'])
            data = data.astype(np.float32) * inunit
            factor = 1
            for base, power in zip(inunit.bases, inunit.powers):
                if (base**power).is_equivalent(u.arcsec**(-2)):
                    conversion =(base**power).to(u.arcsec**(-2)) / base**power
                    data *= conversion
                    factor = u.arcsec**(-2)

            data = data.to(PHOTLAM,
                           equivalencies=u.spectral_density(wave[:, None, None]))

            if factor == 1:    # Normalise to 1 arcsec2 if not a spatial density
                # ..todo: lower needed because "DEG" is not understood, this is ugly
                pixarea = (header['CDELT1'] * u.Unit(header['CUNIT1'].lower()) *
                           header['CDELT2'] * u.Unit(header['CUNIT2'].lower())).to(u.arcsec**2)
                data = data / pixarea.value    # cube is per arcsec2

            data = (data * factor).value

            cube.header['BUNIT'] = 'PHOTLAM/arcsec2'    # ..todo: make this more explicit?

            # The imageplane_utils like to have the spatial WCS in units of "deg". Ensure
            # that the cube is passed on accordingly
            cube.header['CDELT1'] = header['CDELT1'] * u.Unit(header['CUNIT1'].lower()).to(u.deg)
            cube.header['CDELT2'] = header['CDELT2'] * u.Unit(header['CUNIT2'].lower()).to(u.deg)
            cube.header['CUNIT1'] = 'deg'
            cube.header['CUNIT2'] = 'deg'

            # Put on fov wavegrid
            # ..todo: This assumes that we have only one fov. Generalise?
            fov = self.fov_manager.fovs[0]
            wave_min = fov.meta["wave_min"]        # Quantity [um]
            wave_max = fov.meta["wave_max"]

            wave_unit = u.Unit(from_currsys("!SIM.spectral.wave_unit"))
            dwave = from_currsys("!SIM.spectral.spectral_bin_width")  # Not a quantity
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
            cube.header['CTYPE3'] = 'WAVE'
            cube.header['CRPIX3'] = 1
            cube.header['CRVAL3'] = wave_min.value
            cube.header['CDELT3'] = dwave
            cube.header['CUNIT3'] = wave_unit.name

        return source

    def readout(self, filename=None, **kwargs):
        """
        Produces detector readouts for the observed image

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

        hdus = []
        for i, detector_array in enumerate(self.detector_arrays):
            array_effects = self.optics_manager.detector_array_effects
            dtcr_effects = self.optics_manager.detector_effects
            hdu = detector_array.readout(self.image_planes, array_effects,
                                         dtcr_effects, **kwargs)

            if filename is not None and isinstance(filename, str):
                fname = filename
                if len(self.detector_arrays) > 1:
                    fname = str(i) + "_" + filename
                hdu.writeto(fname, overwrite=True)

            hdus += [hdu]

        return hdus

    def set_focus(self, **kwargs):
        self.cmds.update(**kwargs)
        dy = self.cmds.default_yamls
        if len(dy) > 0 and "packages" in dy:
            self.cmds.update(packages=self.default_yamls[0]["packages"])
        rc.__currsys__ = self.cmds


    def shutdown(self):
        '''Shut down the instrument.

        This method closes all open file handles and should be called when the optical train
        is no longer needed.
        '''
        for effect_name in self.effects['name']:
            try:
                self[effect_name]._file.close()
            except AttributeError:
                pass

        self._description = "The instrument has been shut down."


    @property
    def effects(self):
        return self.optics_manager.list_effects()

    def __str__(self):
        return self._description

    def __getitem__(self, item):
        return self.optics_manager[item]

    def __setitem__(self, key, value):
        self.optics_manager[key] = value

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
