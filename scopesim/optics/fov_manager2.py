"""
Influences
Spectral:
    - Red and blue edges of full spectrum
    - Chunks of a large spectral range
Spatial:
    - On-sky borders of Detector Array
    - On-sky borders of Aperture Mask
    - Chunks of large on-sky area
    - Slit mask borders
    - Multiple slits for IFU, mirror-MOS
    - Multiple lenses for fibre-fed MOS

- each Effect should submit a list of volumes

IFU spectroscopy depends on:
    - the tracelist
    - the list of apertures
    - the detector array borders
    - PSF wavelength granularity
    - Atmospheric dispersion

MOS spectroscopy depends on:
    - the tracelist
    - the list of apertures
    - Atmospheric dispersion
    - PSF wavelength granularity
    - the detector array borders

Long slit spectroscopy depends on:
    - the tracelist ra, dec, lam vol --> x, y area
    - the slit aperture
    - the detector array borders
    - PSF wavelength granularity
    - Atmospheric dispersion

Imaging dependent on:
    - Detector array borders
    - PSF wavelength granularity
    - Atmospheric dispersion


"""





from astropy.table import Table
from . import fov_manager_utils as fmu
from ..effects.effects_utils import is_spectroscope
from ..utils import from_currsys


class FOVManager:
    """
    A class to manage the (monochromatic) image windows covering the target

    Parameters
    ----------
    effects : list of Effect objects
        Passed from optics_manager.fov_setup_effects

    kwargs
    ------
    All observation parameters as passed from UserCommands

    """
    def __init__(self, effects=[], **kwargs):
        self.meta = {"area": "!TEL.area",
                     "pixel_scale": "!INST.pixel_scale",
                     "plate_scale": "!INST.plate_scale",
                     "wave_min": "!SIM.spectral.wave_min",
                     "wave_mid": "!SIM.spectral.wave_mid",
                     "wave_max": "!SIM.spectral.wave_max",
                     "chunk_size": "!SIM.computing.chunk_size",
                     "max_segment_size": "!SIM.computing.max_segment_size",
                     "sub_pixel": "!SIM.sub_pixel.flag",
                     "sub_pixel_fraction": "!SIM.sub_pixel.fraction",
                     "preload_fovs": "!SIM.computing.preload_field_of_views"}
        self.meta.update(kwargs)

        self.effects = effects
        if from_currsys(self.meta["preload_fovs"]) is True:
            self._fovs_list = self.generate_fovs_list()
        else:
            self._fovs_list = []

    def generate_fovs_list(self):
        """
        Generates a series of FieldOfViews objects based self.effects

        Returns
        -------
        fovs : list of FieldOfView objects

        """

        self.table = Table(names=["wave_min", "wave_max", 
                                  "sky_left", "sky_right", 
                                  "sky_top", "sky_bottom",
                                  "det_left", "det_right",
                                  "det_top", "det_bottom",
                                  "pixel_scale", "plate_scale",
                                  "sky_rotation", "det_rotation"])

        fovs = []

        return fovs

    @property
    def fovs(self):
        if from_currsys(self.meta["preload_fovs"]) is False:
            self._fovs_list = self.generate_fovs_list()
        return self._fovs_list

    @property
    def fov_footprints(self, which="both"):
        return None

# Spectroscopy FOV setup
# ======================
#
# self.effects --> fov_setup_effects [ApertureList, TraceList, DetectorList]
#
# generate a set of sky headers for each aperture in the list
#
# for each trace in tracelist,
#   note the aperture_id
#   note the image_plane_id
#
#   get the pixel_size [from DetectorList.table["pixsize"]]
#   crawl along trace
#       find the red/blue extreme pixel edge wavelengths
#       get the wavelength of the nearest full pixel jump (x or y)
#
#   interpolate the image plane positions of the trace table for the new
#       wavelength array
#   for each new wavelength
#       get the rotation angle of the fov
#       get the shear of the fov wrt the next fov
#
#       generate an image-plane header for each wavelength
#
# get the combined 3d-shifts
# for each trace
#   interpolate the 3d-shifts for all wavelengths
#   apply the required fov-shift to each fov header


