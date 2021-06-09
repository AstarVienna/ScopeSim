from copy import deepcopy

import numpy as np
from astropy import units as u

from . import image_plane_utils as imp_utils
from .fov import FieldOfView
from .. import effects as efs
from ..effects.effects_utils import get_all_effects
from ..utils import check_keys


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

    Notes
    -----
    Units returned by fov_grid():
    - wavelength: [um]
    - x_shift, y_shift: [deg]

    """
    required_keys = ["wave_min", "wave_mid", "wave_max",
                     "sub_pixel_fraction", "pixel_scale"]
    check_keys(kwargs, required_keys, action="warning")

    effects = get_all_effects(effects, efs.Shift3D)
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
        z_edges  = np.array([z_edges[0]]  + list(z_edges[ii])  + [z_edges[-1]])

    else:
        z_edges = np.array([kwargs["wave_min"], kwargs["wave_max"]])
        x_shifts = np.zeros(2)
        y_shifts = np.zeros(2)

    shift_dict = {"wavelengths": z_edges,
                  "x_shifts": x_shifts / 3600.,   # fov_grid returns [arcsec]
                  "y_shifts": y_shifts / 3600.}   # get_3d_shifts returns [deg]

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
    required_keys = ["wave_min", "wave_max"]
    check_keys(kwargs, required_keys, action="error")

    # get the filter wavelengths first to set (wave_min, wave_max)
    filters = get_all_effects(effects_list, (efs.FilterCurve, efs.FilterWheel))

    wave_bin_edges = [filt.fov_grid(which="waveset", **kwargs)
                      for filt in filters]
    if len(wave_bin_edges) > 0:
        kwargs["wave_min"] = np.max([w[0].value for w in wave_bin_edges])
        kwargs["wave_max"] = np.min([w[1].value for w in wave_bin_edges])
    wave_bin_edges = [[kwargs["wave_min"], kwargs["wave_max"]]]

    if kwargs["wave_min"] > kwargs["wave_max"]:
        raise ValueError("Filter wavelength ranges do not overlap: {}"
                         "".format(wave_bin_edges))

    # ..todo: add in Atmospheric dispersion and ADC here
    for effect_class in [efs.PSF]:
        effects = get_all_effects(effects_list, effect_class)
        for eff in effects:
            waveset = eff.fov_grid(which="waveset", **kwargs)
            if waveset is not None:
                wave_bin_edges += [waveset]

    wave_bin_edges = combine_wavesets(wave_bin_edges)

    if len(wave_bin_edges) == 0:
        wave_bin_edges = [kwargs["wave_min"], kwargs["wave_max"]]

    return wave_bin_edges


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
    # add image plane WCS information for a direct projection

    required_keys = ["pixel_scale", "plate_scale",
                     "chunk_size", "max_segment_size"]
    check_keys(kwargs, required_keys, action="error")

    plate_scale = kwargs["plate_scale"]     # ["/mm]
    pixel_scale = kwargs["pixel_scale"]     # ["/pix]

    # look for apertures
    aperture_effects = get_all_effects(effects, (efs.ApertureMask,
                                                 efs.SlitWheel,
                                                 efs.ApertureList))
    if len(aperture_effects) == 0:
        detector_arrays = get_all_effects(effects, efs.DetectorList)
        if len(detector_arrays) > 0:
            aperture_effects += [detarr.fov_grid(which="edges",
                                                 pixel_scale=pixel_scale)
                                 for detarr in detector_arrays]
        else:
            raise ValueError("No ApertureMask or DetectorList was provided. At "
                             "least one must be passed to make an ImagePlane: "
                             "{}".format(effects))

    # get aperture headers from fov_grid()
    # - for-loop catches mutliple headers from ApertureList.fov_grid()
    hdrs = []
    for ap_eff in aperture_effects:
        # ..todo:: add this functionality to ApertureList effect
        hdr = ap_eff.fov_grid(which="edges", pixel_scale=pixel_scale)
        hdrs += hdr if isinstance(hdr, (list, tuple)) else [hdr]

    # check size of aperture in pixels - split if necessary
    sky_hdrs = []
    for hdr in hdrs:
        if hdr["NAXIS1"] * hdr["NAXIS2"] > kwargs["max_segment_size"]:
            sky_hdrs += imp_utils.split_header(hdr, kwargs["chunk_size"])
        else:
            sky_hdrs += [hdr]

    # ..todo:: Deal with the case that two or more ApertureMasks overlap

    # map the on-sky apertures directly to the image plane using plate_scale
    # - further changes can be made by the individual effects

    # plate_scale ["/mm] --> deg
    # pixel_scale ["/pix]
    # x_sky [deg] / (plate_scale/3600) [deg/mm] --> x_det [mm]

    pixel_size = pixel_scale / plate_scale
    plate_scale_deg = plate_scale / 3600.   # ["/mm] / 3600 = [deg/mm]
    for skyhdr in sky_hdrs:
        x_sky, y_sky = imp_utils.calc_footprint(skyhdr)
        x_det = x_sky / plate_scale_deg
        y_det = y_sky / plate_scale_deg

        dethdr = imp_utils.header_from_list_of_xy(x_det, y_det, pixel_size, "D")
        skyhdr.update(dethdr)

    return sky_hdrs


def get_imaging_fovs(headers, waveset, shifts, **kwargs):
    """
    Returns a list of ``FieldOfView`` objects

    Parameters
    ----------
    headers : list of fits.Header objects
        Headers giving spatial extent of each FOV region

    waveset : list of floats
        [um] N+1 wavelengths for N spectral layers

    shifts : list of tuples
        [deg] x,y shifts w.r.t to the optical axis plane. N shifts for N
        spectral layers

    Returns
    -------
    fovs : list of FieldOfView objects

    """

    shift_waves = shifts["wavelengths"]     # in [um]
    shift_dx = shifts["x_shifts"]           # in [deg]
    shift_dy = shifts["y_shifts"]

    # combine the wavelength bins from 1D spectral effects and 3D shift effects
    if len(shifts["wavelengths"]) > 0:
        mask = (shift_waves > np.min(waveset)) * (shift_waves < np.max(waveset))
        waveset = combine_wavesets([waveset, shift_waves[mask]])

    counter = 0
    fovs = []

    print("Preparing {} FieldOfViews".format((len(waveset)-1)*len(headers)),
          flush=True)

    for ii in range(len(waveset) - 1):
        for hdr in headers:
            # add any pre-instrument shifts to the FOV sky coords
            wave_mid = 0.5 * (waveset[ii] + waveset[ii+1])
            x_shift = np.interp(wave_mid, shift_waves, shift_dx)
            y_shift = np.interp(wave_mid, shift_waves, shift_dy)

            fov_hdr = deepcopy(hdr)
            fov_hdr["CRVAL1"] += x_shift      # headers are in [deg]
            fov_hdr["CRVAL2"] += y_shift

            # define the wavelength range for the FOV
            waverange = [waveset[ii], waveset[ii + 1]]

            # Make the FOV
            fov = FieldOfView(fov_hdr, waverange, id=counter, **kwargs)
            fovs += [fov]
            counter += 1

    return fovs


def get_spectroscopy_headers(effects, **kwargs):

    required_keys = ["pixel_scale", "plate_scale",
                     "wave_min", "wave_max"]
    check_keys(kwargs, required_keys, action="error")

    surface_list_effects = get_all_effects(effects, (efs.SurfaceList,
                                                     efs.FilterWheel,
                                                     efs.FilterCurve))
    detector_list_effects = get_all_effects(effects, efs.DetectorList)
    spec_trace_effects = get_all_effects(effects, efs.SpectralTraceList)
    aperture_effects = get_all_effects(effects, (efs.ApertureList,
                                                 efs.SlitWheel,
                                                 efs.ApertureMask))

    if len(surface_list_effects) > 0:
        waves = surface_list_effects[0].fov_grid(which="waveset")
        if len(waves) == 2:
            kwargs["wave_min"] = np.max([waves[0].value, kwargs["wave_min"]])
            kwargs["wave_max"] = np.min([waves[1].value, kwargs["wave_max"]])

    if len(detector_list_effects) > 0:
        implane_hdr = detector_list_effects[0].image_plane_header
    elif len(spec_trace_effects) > 0:
        implane_hdr = spec_trace_effects[0].image_plane_header
    else:
        raise ValueError("Missing a way to determine the image plane size")

    # ..todo: deal with multiple trace lists
    if len(spec_trace_effects) != 1:
        raise ValueError("More than one SpectralTraceList was found: {}"
                         "".format(spec_trace_effects))
    spec_trace = spec_trace_effects[0]

    sky_hdrs = []
    for ap_eff in aperture_effects:
        # if ApertureList, a list of ApertureMask headers is returned
        # If ApertureMask, a single header is returned
        hdr = ap_eff.fov_grid(which="edges", pixel_scale=kwargs["pixel_scale"])
        sky_hdrs += hdr if isinstance(hdr, (list, tuple)) else [hdr]

    fov_headers = [spec_trace.fov_grid(which="edges",
                                       sky_header=sky_hdr,
                                       det_header=implane_hdr,
                                       wave_min=kwargs["wave_min"],
                                       wave_max=kwargs["wave_max"],
                                       pixel_scale=kwargs["pixel_scale"],
                                       plate_scale=kwargs["plate_scale"]
                                       )
                   for sky_hdr in sky_hdrs]
    fov_headers = [hdr for hdr_list in fov_headers for hdr in hdr_list]
    # ..todo: check that each header is not larger than chunk_size

    return fov_headers


def get_spectroscopy_fovs(headers, shifts, effects=[], **kwargs):

    shift_waves = shifts["wavelengths"]     # in [um]
    shift_dx = shifts["x_shifts"]           # in [deg]
    shift_dy = shifts["y_shifts"]

    print("Preparing {} FieldOfViews".format(len(headers)), flush=True)

    apertures = get_all_effects(effects, (efs.ApertureList, efs.ApertureMask))
    masks = [ap.fov_grid(which="masks") for ap in apertures]
    mask_dict = {}
    for mask in masks:
        if isinstance(mask, dict):
            mask_dict.update(mask)
        elif isinstance(mask, np.ndarray):
            mask_dict[len(mask_dict)] = mask

    fovs = []
    for ii, hdr in enumerate(headers):
        # add any pre-instrument shifts to the FOV sky coords
        wave_mid = hdr["WAVE_MID"]
        x_shift = np.interp(wave_mid, shift_waves, shift_dx)
        y_shift = np.interp(wave_mid, shift_waves, shift_dy)

        fov_hdr = deepcopy(hdr)
        fov_hdr["CRVAL1"] += x_shift      # headers are in [deg]
        fov_hdr["CRVAL2"] += y_shift

        # Make the FOV
        fov = FieldOfView(fov_hdr, waverange=[hdr["WAVE_MIN"], hdr["WAVE_MAX"]],
                          **kwargs)
        fov.meta["distortion"]["rotation"] = hdr["ROTANGD"]
        fov.meta["distortion"]["shear"] = hdr["SKEWANGD"]
        fov.meta["conserve_image"] = hdr["IMG_CONS"]
        fov.meta["fov_id"] = ii
        fov.meta["aperture_id"] = hdr["APERTURE"]

        # .. todo: get these masks working
        # there needs to be fov_grid(which="mask") in ApertureList/Mask
        # fov.mask = mask_dict[hdr["APERTURE"]]
        fovs += [fov]

    return fovs


def combine_wavesets(waveset_list):
    """
    Joins and sorts several sets of wavelengths into a single 1D array

    Parameters
    ----------
    waveset_list : list
        A group of wavelength arrays or lists

    Returns
    -------
    wave_set : np.ndarray
        Combined set of wavelengths

    """
    wave_set = []
    for wbe in waveset_list:        # wbe = waveset bin edges
        wbe = wbe.value if isinstance(wbe, u.Quantity) else wbe
        wave_set += list(wbe)
    # ..todo:: set variable in !SIM.computing for rounding to the 7th decimal
    wave_set = np.unique(np.round(np.sort(wave_set, kind="stable"), 7))

    return wave_set
