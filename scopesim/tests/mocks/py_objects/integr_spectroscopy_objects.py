import numpy as np
from astropy.io import fits
from astropy.io.fits import table_to_hdu
from astropy.table import Table
from astropy import units as u
from synphot import SourceSpectrum
from synphot.units import PHOTLAM
from synphot.models import Empirical1D

from scopesim import effects as efs
from scopesim.source.source import Source


def mock_aperture_list():
    """
    Two rectangular slits (dimensions 10"x1"), above and below the optical axis
    """
    array_dict = {"id": [0, 1], "left": [-5, -5], "right": [5, 5],
                  "top": [1.5, -0.5], "bottom": [0.5, -1.5],
                  "angle": [0, 0], "conserve_image": [True, True],
                  "shape": ["rect", "rect"]}
    kwargs = {"x_unit": "arcsec", "y_unit": "arcsec", "angle_unit": "deg",
              "pixel_scale": 0.1}
    ap_list = efs.ApertureList(array_dict=array_dict, **kwargs)

    return ap_list


def mock_aperture_list_mixed():
    """
    Two rectangular slits (dimensions 10"x1"), above and below the optical axis
    """
    array_dict = {"id": [0, 1], "left": [-0.5, -5], "right": [0.5, 5],
                  "top": [1.5, -0.5], "bottom": [0.5, -1.5],
                  "angle": [0, 0], "conserve_image": [False, True],
                  "shape": ["round", "rect"]}
    kwargs = {"x_unit": "arcsec", "y_unit": "arcsec", "angle_unit": "deg",
              "pixel_scale": 0.1, "no_mask": False}
    ap_list = efs.ApertureList(array_dict=array_dict, **kwargs)

    return ap_list


def mock_aperture_list_single():
    """
    Two rectangular slits (dimensions 10"x1"), above and below the optical axis
    """
    array_dict = {"id": [0], "left": [-5], "right": [5],
                  "top": [-0.5], "bottom": [-1.5],
                  "angle": [0], "conserve_image": [True],
                  "shape": ["rect"]}
    kwargs = {"x_unit": "arcsec", "y_unit": "arcsec", "angle_unit": "deg",
              "pixel_scale": 0.1}
    ap_list = efs.ApertureList(array_dict=array_dict, **kwargs)

    return ap_list


def mock_spectral_trace_list_single():
    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header.update({"ECAT": 1, "EDATA": 2})
    idx_hdu = table_to_hdu(Table(names=['description', 'extension_id',
                                        'aperture_id', 'image_plane_id'],
                                 data=[["0"], [2], [0], [0]]))
    n = 11
    names = ["wavelength", "s0", "s1", "x0", "x1", "y0", "y1"]
    wA, wB = np.linspace(1., 2., n), np.linspace(1.8, 2.4, n)
    s0, s1 = -5 * np.ones(n), 5 * np.ones(n)
    x0 = np.linspace(-25, 25, n)
    x1 = x0 + 10
    y0 = np.linspace(-50, 50, n)
    y1 = y0 +2

    trace_hdu = table_to_hdu(Table(names=names, data=[wA, s0, s1, x0, x1, y0, y1]))
    hdu_list = fits.HDUList([pri_hdu, idx_hdu, trace_hdu])
    spt_list = efs.SpectralTraceList(hdulist=hdu_list)

    return spt_list


def mock_spectral_trace_list():
    """
    2 apertures, each with two traces. The main traces [0_A, 1_A] cover the
    waverange [1, 2]um, and are on the left of the image plane, spanning two
    detector. The secondary traces cover [1.8, 2.4]um and are right of centre
    """

    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header.update({"ECAT": 1, "EDATA": 2})
    idx_hdu = table_to_hdu(Table(names=['description', 'extension_id',
                                        'aperture_id', 'image_plane_id'],
                                 data=[["0_A", "0_B", "1_A", "1_B"],
                                       [2, 3, 4, 5], [0, 0, 1, 1], [0]*4]))
    names = ["wavelength", "s0", "s1", "x0", "x1", "y0", "y1"]
    wA, wB = [1.0, 2.0], [1.8, 2.4]
    s0, s1 = [-5, -5], [5, 5]
    trace_hdus = [table_to_hdu(Table(names=names,               # Trace 0_A
                                     data=[wA, s0, s1,
                                           [-50, -50], [-40, -40],
                                           [-50, 50], [-50, 50]])),
                  table_to_hdu(Table(names=names,               # Trace 0_B
                                     data=[wB, s0, s1,
                                           [10, 10], [20, 20],
                                           [-20, 20], [-20, 20]])),
                  table_to_hdu(Table(names=names,               # Trace 1_A
                                     data=[wA, s0, s1,
                                           [-30, -30], [-20, -20],
                                           [-50, 50], [-50, 50]])),
                  table_to_hdu(Table(names=names,               # Trace 1_B
                                     data=[wB, s0, s1,
                                           [30, 30], [40, 40],
                                           [-20, 20], [-20, 20]]))]
    hdu_list = fits.HDUList([pri_hdu, idx_hdu] + trace_hdus)
    kwargs = {}

    spt_list = efs.SpectralTraceList(hdulist=hdu_list, **kwargs)

    return spt_list


def mock_spectral_trace_list_shear():
    """
    2 apertures, each with two traces. The main traces [0_A, 1_A] cover the
    waverange [1, 2]um, and are on the left of the image plane, spanning two
    detectors. The secondary traces cover [1.8, 2.4]um and are right of centre
    """

    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header.update({"ECAT": 1, "EDATA": 2})
    idx_hdu = table_to_hdu(Table(names=['description', 'extension_id',
                                        'aperture_id', 'image_plane_id'],
                                 data=[["0_A", "0_B", "1_A", "1_B"],
                                       [2, 3, 4, 5], [1, 1, 1, 1], [0]*4]))
    n = 11
    names = ["wavelength", "s0", "s1", "x0", "x1", "y0", "y1"]
    wA, wB = np.linspace(1., 2., n), np.linspace(1.8, 2.4, n)
    s0, s1 = -5 * np.ones(n), 5 * np.ones(n)

    # diagonal
    y_A = np.linspace(-50, 50, n)
    x0_0A = np.linspace(-50, -40, n)
    x1_0A = x0_0A + 10

    # straight
    x0_1A = -20 * np.ones(n)
    x1_1A = x0_1A + 10

    # sheared
    x0_0B = np.linspace(0, 10, n)
    y0_0B = np.linspace(-20, 20, n)
    x1_0B = x0_0B + 10
    y1_0B = y0_0B + 10

    # horizontal
    x0_1B = np.linspace(20, 50, n)
    y0_1B = -5 * np.ones(n)
    x1_1B = x0_1B
    y1_1B = y0_1B + 10

    trace_hdus = [table_to_hdu(Table(names=names,               # Trace 0_A
                                     data=[wA, s0, s1, x0_0A, x1_0A, y_A, y_A])),
                  table_to_hdu(Table(names=names,               # Trace 0_B
                                     data=[wB, s0, s1, x0_0B, x1_0B, y0_0B, y1_0B])),
                  table_to_hdu(Table(names=names,               # Trace 1_A
                                     data=[wA, s0, s1, x0_1A, x1_1A, y_A, y_A])),
                  table_to_hdu(Table(names=names,               # Trace 1_B
                                     data=[wB, s0, s1, x0_1B, x1_1B, y0_1B, y1_1B]))]
    hdu_list = fits.HDUList([pri_hdu, idx_hdu] + trace_hdus)
    kwargs = {}

    spt_list = efs.SpectralTraceList(hdulist=hdu_list, **kwargs)

    return spt_list


def mock_detector_list():
    """
    Four 512x512 detectors around the optical axis, seperated by 100px.
    Pixel size is 0.1mm
    """
    array_dict = {"id": [0, 1, 2, 3], "pixsize": [0.1]*4, "angle": [0.]*4,
                  "gain": [1.]*4, "xhw": [25.6]*4, "yhw": [25.6]*4,
                  "x_cen": [-30.6, -30.6, 30.6, 30.6],
                  "y_cen": [30.6, -30.6, -30.6, 30.6]}
    kwargs = {"x_cen_unit": "mm", "y_cen_unit": "mm",
              "xhw_unit": "mm", "yhw_unit": "mm", "pixsize_unit": "mm",
              "angle_unit": "deg", "gain_unit": "electron/adu",
              "image_plane_id": 0}
    det_list = efs.DetectorList(array_dict=array_dict, **kwargs)

    return det_list


def mock_config_yaml():
    config = {"pixel_scale": 0.1,       # arcsec / pix
              "plate_scale": 1,         # arcsec / mm
              "wave_min": 1.0,
              "wave_mid": 1.6,
              "wave_max": 2.5}
    return config


def mock_point_source_object():
    """
    Three point sources, 1 up and 2 below the optical axis, fit inside the
    above defined apertures
    """

    wave = [1.0, 2.5] * u.um
    flat = SourceSpectrum(Empirical1D, points=wave,
                          lookup_table=[100, 100] * PHOTLAM)
    down = SourceSpectrum(Empirical1D, points=wave,
                          lookup_table=[60, 1] * PHOTLAM)
    up   = SourceSpectrum(Empirical1D, points=wave,
                          lookup_table=[1, 60] * PHOTLAM)

    src = Source(x=[0, -2, 2], y=[1, -1, -1], ref=[0, 1, 2],
                 spectra=[flat, up, down])

    return src


def mock_extended_source_object():
    from scipy.misc import face
    im = face()[::-1, :, 0]
    im = im / np.max(im)
    hdr = {"CDELT1": 0.1, "CDELT2": 0.1,
           "CRVAL1": 0.,  "CRVAL2": 0.,
           "CRPIX1": im.shape[1]/2, "CRPIX2": im.shape[0]/2,
           "CTYPE1": "LINEAR", "CTYPE2": "LINEAR",
           "CUNIT1": "ARCSEC", "CUNIT2": "ARCSEC"}
    hdu = fits.ImageHDU(data=im, header=fits.Header(hdr))
    spec = SourceSpectrum(Empirical1D, points=[1.0, 2.5] * u.um,
                          lookup_table=[1, 1] * PHOTLAM)
    src = Source(image_hdu=hdu, spectra=spec)

    return src


def mock_gauss_psf():
    return efs.GaussianDiffractionPSF(diameter=1*u.m, convolve_mode="same")


def mock_3d_shift():
    array_dict = {"wavelength": [1.0, 2.5], "dx": [0, 0], "dy": [-1, 1]}
    return efs.Shift3D(array_dict=array_dict)
