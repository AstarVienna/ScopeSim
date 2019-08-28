from copy import deepcopy

from .. import rc
from ..commands.user_commands import UserCommands
from .optics_manager import OpticsManager
from .fov_manager import FOVManager
from .image_plane import ImagePlane
from ..detector import DetectorArray


class OpticalTrain:
    """
    The main class for controlling a simulation

    Parameters
    ----------
    cmds : UserCommands

    """

    def __init__(self, cmds=None):

        self.cmds = cmds
        self.optics_manager = None
        self.fov_manager = None
        self.image_plane = None
        self.detector_array = None
        self.yaml_dicts = None

        if cmds is not None:
            self.load(cmds)

    def load(self, user_commands):
        """
        (Re)Loads an OpticalTrain with a new set of UserCommands

        Parameters
        ----------
        user_commands : UserCommands

        """

        if not isinstance(user_commands, UserCommands):
            raise ValueError("user_commands must be a UserCommands object: "
                             "{}".format(type(user_commands)))

        rc.__currsys__ = user_commands
        self.yaml_dicts = rc.__currsys__.yaml_dicts
        self.optics_manager = OpticsManager(self.yaml_dicts)
        self.update()

    def update(self, **kwargs):
        """
        Update the user-defined parameters and remake the main internal classes

        Parameters
        ----------
        kwargs : **dict
            Any keyword-value pairs from a config file

        """
        self.optics_manager.update(**kwargs)
        opt_man = self.optics_manager
        self.fov_manager = FOVManager(opt_man.fov_setup_effects, **kwargs)
        self.image_plane = ImagePlane(opt_man.image_plane_header, **kwargs)
        self.detector_array = DetectorArray(opt_man.detector_setup_effects,
                                            **kwargs)

    def observe(self, orig_source, **kwargs):
        """
        Main controlling method for observing ``Source`` objects

        Parameters
        ----------
        orig_source : Source
        kwargs : **dict
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

        """

        self.update(**kwargs)

        source = deepcopy(orig_source)

        # [1D - transmisison curves]
        for effect in self.optics_manager.source_effects:
            source = effect.apply_to(source)

        # [3D - Atmospheric shifts, PSF, NCPAs, Grating shift/distortion]
        for fov in self.fov_manager.fovs:
            fov.extract_from(source)
            for effect in self.optics_manager.fov_effects:
                fov = effect.apply_to(fov)

            if fov.hdu.data is None:
                fov.view()
            self.image_plane.add(fov.hdu, wcs_suffix="D")

        # [2D - Vibration, flat fielding, chopping+nodding]
        for effect in self.optics_manager.image_plane_effects:
            self.image_plane = effect.apply_to(self.image_plane)

    def readout(self, filename=None, **kwargs):
        """

        Parameters
        ----------
        filename : str
        kwargs

        Returns
        -------
        hdu : fits.HDUList

        Notes
        -----
        - Apply detector plane (0D, 2D) effects - z_order = 500..599

        """

        hdu = None
        if self.detector_array is not None:
            dtcr_effects = self.optics_manager.detector_effects
            hdu = self.detector_array.readout(self.image_plane, dtcr_effects,
                                              **kwargs)

        if filename is not None and isinstance(filename, str):
            hdu.writeto(filename, overwrite=True)

        return hdu
