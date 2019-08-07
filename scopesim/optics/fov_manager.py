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

import numpy as np
from astropy import units as u

import scopesim.effects.apertures
from .image_plane_utils import header_from_list_of_xy
from .fov import FieldOfView
from .. import effects as efs
from ..effects import Shift3D
from ..effects import FilterCurve, PSF, SpectralTraceList
from ..effects.effects_utils import get_all_effects, is_spectroscope
from ..utils import check_keys, from_currsys, quantify


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
        self.meta = {"pixel_scale": "!INST.pixel_scale",
                     "wave_min": "!SIM.spectral.wave_min",
                     "wave_mid": "!SIM.spectral.wave_mid",
                     "wave_max": "!SIM.spectral.wave_max",
                     "sub_pixel_fraction": "!SIM.sub_pixel.fraction",
                     "chunk_size": "!SIM.computing.chunk_size"}
        self.meta.update(kwargs)

        self.effects = effects
        self._fovs_list = []

    def generate_fovs_list(self):
        """
        Generates a series of FieldOfViews objects based self.effects

        Returns
        -------
        fovs : list of FieldOfView objects

        """
        self.meta = from_currsys(self.meta)

        if is_spectroscope(self.effects):
            shifts  = get_3d_shifts(self.effects, **self.meta)
            waveset = get_spectroscopy_waveset(self.effects, **self.meta)
            headers = get_spectroscopy_headers(self.effects, **self.meta)
            fovs    = get_spectroscopy_fovs(headers, waveset, shifts)
        else:
            shifts  = get_3d_shifts(self.effects, **self.meta)
            waveset = get_imaging_waveset(self.effects, **self.meta)
            headers = get_imaging_headers(self.effects, **self.meta)
            fovs    = get_imaging_fovs(headers, waveset, shifts)

        return fovs

    @property
    def fovs(self):
        self._fovs_list = self.generate_fovs_list()
        return self._fovs_list

    @property
    def fov_footprints(self):
        return None


def get_3d_shifts(effects, **kwargs):
    """
    Returns the total 3D shifts (x,y,lam) from a series of Shift3D objects

    Parameters
    ----------
    effects : list of Shift3D effects

    Returns
    -------
    shift_dict : dict
        returns the x, y shifts for each wavelength in the fov_grid, where
        fov_grid contains the edge wavelengths for each spectral layer

    """
    required_keys = ["wave_min", "wave_mid", "wave_max",
                     "sub_pixel_fraction", "pixel_scale"]
    check_keys(kwargs, required_keys, action="warning")

    effects = get_all_effects(effects, Shift3D)
    if len(effects) > 0:
        # shifts = [[waves], [x_shifts], [y_shifts]]
        shifts = [eff.fov_grid(which="shifts", **kwargs) for eff in effects]

        old_bin_edges = [shift[0] for shift in shifts if len(shift[0]) >= 2]
        new_bin_edges = np.unique(np.sort(np.concatenate(old_bin_edges),
                                          kind="stable"))

        x_shifts = np.zeros(len(new_bin_edges))
        y_shifts = np.zeros(len(new_bin_edges))
        # .. todo:: replace the 1e-7 with a variable in !SIM
        for shift in shifts:
            if not np.all(np.abs(shift[1]) < 1e-7):
                x_shifts += np.interp(new_bin_edges, shift[0], shift[1])
            if not np.all(np.abs(shift[2]) < 1e-7):
                y_shifts += np.interp(new_bin_edges, shift[0], shift[2])

        # After adding all the shifts together, work out a new wavelength set
        z_edges = np.copy(new_bin_edges)
        z_shifts = (x_shifts**2 + y_shifts**2)**0.5        # in arcsec
        step_size = kwargs["pixel_scale"] * kwargs["sub_pixel_fraction"]
        z_steps = (z_shifts / step_size).astype(int)
        # find where the shift is larger than a sub pixel fraction size
        ii = np.where(np.diff(z_steps) != 0)[0]

        x_shifts = np.array([x_shifts[0]] + list(x_shifts[ii]) + [x_shifts[-1]])
        y_shifts = np.array([y_shifts[0]] + list(y_shifts[ii]) + [y_shifts[-1]])
        z_edges = np.array([z_edges[0]] + list(z_edges[ii]) + [z_edges[-1]])

    else:
        z_edges = np.array([kwargs["wave_min"], kwargs["wave_max"]])
        x_shifts = np.zeros(2)
        y_shifts = np.zeros(2)

    shift_dict = {"wavelengths": z_edges,
                  "x_shifts": x_shifts,
                  "y_shifts": y_shifts}

    return shift_dict


def get_imaging_waveset(effects_list, **kwargs):
    """
    Returns the edge wavelengths for the spectral layers needed for simulation

    Parameters
    ----------
    effects_list : list of Effect objects

    Returns
    -------
    wave_bin_edges : list
        [um] list of wavelengths

    """
    # get the filter wavelengths first to set (wave_min, wave_max)
    filters = get_all_effects(effects_list, FilterCurve)

    wave_bin_edges = [filt.fov_grid(which="waveset", **kwargs)
                      for filt in filters]
    if len(wave_bin_edges) > 0:
        kwargs["wave_min"] = np.min(wave_bin_edges)
        kwargs["wave_max"] = np.max(wave_bin_edges)

    for effect_class in (PSF, SpectralTraceList):
        effects = get_all_effects(effects_list, effect_class)
        for eff in effects:
            wave_bin_edges += [eff.fov_grid(which="waveset", **kwargs)]

    wave_set = []
    for wbe in wave_bin_edges:
        wbe = wbe.value if isinstance(wbe, u.Quantity) else wbe
        wave_set += list(wbe)
    # ..todo:: set variable in !SIM for the rounding to the 7th decimal
    wave_set = np.unique(np.round(np.sort(wave_set, kind="stable"), 7))

    return wave_set


def get_imaging_headers(effects, **kwargs):
    """
    Returns a list of Header objects for each of the FieldOfVIew objects

    Parameters
    ----------
    effects : list of Effect objects
        Should contain all effects which return a header defining the extent
        of their spatial coverage

    Returns
    -------
    hdrs : list of Header objects

    Notes
    -----
    FOV headers use the return values from the ``<Effect>.fov_grid()``
    method. The ``fov_grid`` dict must contain the entry ``edges``

    This may change in future versions of ScopeSim

    """
    # look for apertures,
    #   if no apertures, look for detectors, make aperture from detector
    # get aperture headers via fov_grid()
    # check size of aperture in pixels
    #   if larger than max_chunk_size, split into smaller headers

    required_keys = ["pixel_scale", "chunk_size", "max_segment_size"]
    check_keys(kwargs, required_keys, action="warning")

    # look for apertures
    aperture_effects = get_all_effects(effects, (efs.ApertureMask,
                                                 efs.ApertureList))
    if len(aperture_effects) == 0:
        detector_arrays = get_all_effects(effects, efs.DetectorList)
        if len(detector_arrays) > 0:
            # ..todo:: add this functionality to DetectorList effect
            pixel_scale = kwargs["pixel_scale"]
            aperture_effects += [detarr.fov_grid(which="edges",
                                                 pixel_scale=pixel_scale)
                                 for detarr in detector_arrays]
        else:
            raise ValueError("No ApertureMask or DetectorList was provided. At "
                             "least one must be passed to make an ImagePlane: "
                             "{}".format(effects))

    # get aperture headers from fov_grid()
    hdrs = []
    for ap_eff in aperture_effects:
        # - for loop catches mutliple headers from ApertureList.fov_grid()
        # ..todo:: add this functionality to ApertureList effect
        hdr = ap_eff.fov_grid(which="edges")
        hdrs += list(hdr) if isinstance(hdr, (list, tuple)) else [hdr]

    # check size of aperture in pixels - split if necessary
    new_hdrs = []
    for hdr in hdrs:
        if hdr["NAXIS1"] * hdr["NAXIS2"] > kwargs["max_segment_size"]:
            new_hdrs += split_header(hdr, kwargs["chunk_size"])
        else:
            new_hdrs += [hdr]

    return new_hdrs


def make_aperture_from_detector(detector_arrays, pixel_scale):
    pass


def split_header(hdr, chunk_size):
    pass


    #
    # pixel_scale_deg = kwargs["pixel_scale"] / 3600.         # " -> deg
    # pixel_size = detector_arrays[0].image_plane_header["CDELT1D"]  # mm
    # deg2mm = pixel_size / pixel_scale_deg
    #
    # sky_edges = []
    # if len(aperture_masks) > 0:
    #     sky_edges += [apm.fov_grid(which="edges") for apm in aperture_masks]
    # elif len(detector_arrays) > 0:
    #     edges = detector_arrays[0].fov_grid(which="edges",
    #                                         pixel_scale=pixel_scale_deg)
    #     sky_edges += [edges]
    # else:
    #     raise ValueError("No ApertureMask or DetectorList was provided. At "
    #                      "least a DetectorList object must be passed: {}"
    #                      "".format(effects))
    #
    # width = kwargs["chunk_size"] * pixel_scale_deg
    # hdrs = []
    # for xy_sky in sky_edges:
    #     x0, y0 = min(xy_sky[0]), min(xy_sky[1])
    #     x1, y1 = max(xy_sky[0]), max(xy_sky[1])
    #     for xi in np.arange(x0, x1, width):
    #         for yi in np.arange(y0, y1, width):
    #             # xii = np.array([xi, xi + min(width, x1-xi)])
    #             # yii = np.array([yi, yi + min(width, y1-yi)])
    #             xii = np.array([xi, xi + width])
    #             yii = np.array([yi, yi + width])
    #             hdr_sky = header_from_list_of_xy(xii, yii, pixel_scale_deg)
    #             hdr_mm  = header_from_list_of_xy(xii * deg2mm, yii * deg2mm,
    #                                              pixel_size, "D")
    #             hdr_mm.update(hdr_sky)
    #             hdrs += [hdr_mm]
    #
    # return hdrs


def get_imaging_fovs(headers, waveset, shifts):
    """
    Returns a list of ``FieldOfView`` objects

    Parameters
    ----------
    headers : list of Header objects

    waveset : list of floats
        [um] N+1 wavelengths for N spectral layers

    shifts : list of tuples
        [deg] x,y shifts w.r.t to the optical axis plane. N shifts for N
        spectral layers

    Returns
    -------
    fovs : list of FieldOfView objects

    """

    if len(shifts["wavelengths"]) > len(waveset):
        waveset = shifts["wavelengths"]

    # ..todo: add the shifts in somehow

    counter = 0
    fovs = []
    for ii in range(len(waveset) - 1):
        for hdr in headers:
            waverange = [waveset[ii], waveset[ii + 1]]
            fov = FieldOfView(hdr, waverange)
            fov.meta["id"] = counter
            fovs += [fov]

    return fovs


def get_spectroscopy_waveset(effects, **kwargs):
    # TraceList
    # for each trace dlam along the trace centre in increments
    #   of SIM_SUB_PIXEL_FRACTION
    # Must be accompanied by an ApertureList
    pass


def get_spectroscopy_headers(effects, **kwargs):
    # ApertureList
    # For each Trace set the sky header to the aperture footprint
    #   plus any shifts from AtmosphericDispersion
    # Set the Image plane footprint centred on the image plane
    #   position
    pass


def get_spectroscopy_fovs(fields, waveset, shifts):
    pass
