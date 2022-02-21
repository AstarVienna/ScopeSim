# 1. Find the Wavelength range
# Build from edges of throughput curve

# 2. Find the wavelength bins
# If TraceList and Aperture list, then Spectroscopy
# TraceList
# for each trace dlam along the trace centre in increments
#   of SIM_SUB_PIXEL_FRACTION
# Must be accompanied by an ApertureList

# If not, then imaging
# PSF core increase (atmo, ncpas)
# If from a files, what is the bin size?
# If analytic, dlam between a FWHM or SIM_SUB_PIXEL_FRACTION
# ADC + AD shifts
# dlam between shift of SIM_SUB_PIXEL_FRACTION

# 3. Find the spatial range
# If Spectroscopy
# ApertureList
# For each Trace set the sky header to the aperture footprint
#   plus any shifts from AtmosphericDispersion
# Set the Image plane footprint centred on the image plane
#   position

# If Imaging
# DetectorList, or ApertureMask, plus any shift from
#   AtmosphericDispersion

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

        self.meta = from_currsys(self.meta)
        fov_meta_keys = ["area", "sub_pixel"]
        fov_meta = {key: self.meta[key] for key in fov_meta_keys}

        if is_spectroscope(self.effects):
            headers = fmu.get_spectroscopy_headers(self.effects, **self.meta)
            shifts = fmu.get_3d_shifts(self.effects, **self.meta)
            fovs = fmu.get_spectroscopy_fovs(headers, shifts, self.effects,
                                             **fov_meta)
        else:
            headers = fmu.get_imaging_headers(self.effects, **self.meta)
            waveset = fmu.get_imaging_waveset(self.effects, **self.meta)
            shifts = fmu.get_3d_shifts(self.effects, **self.meta)
            fovs = fmu.get_imaging_fovs(headers, waveset, shifts, **fov_meta)

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


