# actually for Source
import pytest
from pytest import approx

import os
from copy import deepcopy

import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table
from astropy import units as u
from astropy import wcs

from synphot import SourceSpectrum, SpectralElement
from synphot.models import Empirical1D
from synphot.units import PHOTLAM

import scopesim as sim
from scopesim.source import source_utils
from scopesim.source.source import Source

from scopesim.optics.image_plane import ImagePlane
from scopesim.utils import convert_table_comments_to_dict

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


MOCK_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        "../mocks/files/"))
sim.rc.__search_path__.insert(0, MOCK_DIR)

PLOTS = True


def make_dummy_cube(scale, wave_unit, ref_wave, wave_step, wave_type, bunit):

    if isinstance(scale, u.Quantity) is False:
        scale = scale * u.arcsec
    if isinstance(wave_unit, u.core.Unit) is False:
        wave_unit = u.AA
    if isinstance(ref_wave, u.Quantity) is False:
        ref_wave = ref_wave * u.AA
    if isinstance(wave_step, u.Quantity) is False:
        wave_step = wave_step * u.AA

    data = np.ones(shape=(1000, 200, 200))
    header = fits.Header(dict(NAXIS=3, NAXIS1=data.shape[0] + 1,
                              NAXIS2=data.shape[1] + 1,
                              NAXIS3=data.shape[2] + 1,
                              CRPIX1=data.shape[0] // 2,
                              CRPIX2=data.shape[0] // 2,
                              CRPIX3=1,
                              CRVAL1=0,
                              CRVAL2=0,
                              CRVAL3=ref_wave.to(wave_unit).value,
                              CDELT1=-1 * scale.to(u.deg).value,
                              CDELT2=scale.to(u.deg).value,
                              CDELT3=wave_step.to(wave_unit).value,
                              CUNIT1="DEG",
                              CUNIT2="DEG",
                              CUNIT3=wave_unit.to_string(),
                              CTYPE1='RA---TAN',
                              CTYPE2='DEC--TAN',
                              CTYPE3=wave_type,
                              BUNIT=bunit))

    hdu = fits.PrimaryHDU(data=data, header=header)

    return hdu


def test_source_cube_angstrom():
    """
    Testing basics
    """

    cube = make_dummy_cube(scale=0.2*u.arcsec,
                           wave_unit=u.AA, ref_wave=1000*u.AA, wave_step=1*u.AA, wave_type="WAVE",
                           bunit="erg / (s angstrom cm2 arcsec2)")

    src_cube = Source(cube=cube)

    if PLOTS:
        src_cube.plot()
        plt.show()

    assert isinstance(src_cube.fields[0], fits.PrimaryHDU)


def test_source_hdu_list():
    """
    Testing basics
    """

    cube = make_dummy_cube(scale=0.2*u.arcsec,
                           wave_unit=u.AA, ref_wave=1000*u.AA, wave_step=1*u.AA, wave_type="WAVE",
                           bunit="erg / (s angstrom cm2 arcsec2)")

    cube_list = fits.HDUList(hdus=[cube, cube])

    src_cube = Source(cube=cube_list, ext=1)

    if PLOTS:
        src_cube.plot()
        plt.show()

    assert isinstance(src_cube.fields[0], fits.PrimaryHDU)


def test_source_from_file():
    """
    Testing basics
    """

    cube = make_dummy_cube(scale=0.2*u.arcsec,
                           wave_unit=u.AA, ref_wave=1000*u.AA, wave_step=1*u.AA, wave_type="WAVE",
                           bunit="erg / (s angstrom cm2 arcsec2)")

    filename = os.path.join(MOCK_DIR, "cube.fits")
    cube.writeto(filename)

    src_cube = Source(cube=filename)

    if PLOTS:
        src_cube.plot()
        plt.show()

    assert isinstance(src_cube.fields[0], fits.PrimaryHDU)

