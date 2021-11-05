"""
Warning: To test the creation of plots, if PLOTS=False, all tests will pass regardless
"""

import os
import pytest
from pytest import approx

import numpy as np
from astropy import units as u
from astropy.io import fits

from scopesim import rc
from scopesim.effects import psfs, electronic, shifts

from scopesim.tests.mocks.py_objects import fov_objects as fovobj
from scopesim.tests.mocks.py_objects.fov_objects import _centre_fov

from scopesim.tests.mocks.py_objects import source_objects as srcobj

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False  # if FALSE all tests will pass regardless

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


class TestPsfEffects:
    def test_VibrationPsf_plot(self):
        psf = psfs.Vibration(fwhm=0.03, pixel_scale=0.005)
        fov = fovobj._centre_micado_fov()
        if PLOTS:
            psf.plot(fov)
            plt.show()

        assert True

    def test_NonCommonPathAberrationPsf_plot(self):
        psf = psfs.NonCommonPathAberration(pixel_scale=0.05)
        fov = fovobj._centre_micado_fov()

        if PLOTS:
            psf.plot(fov)
            plt.show()

        assert True

    def test_SeeingPSF_plot(self):
        psf = psfs.SeeingPSF(fwhm=0.8)
        fov = fovobj._centre_micado_fov()
        if PLOTS:
            psf.plot(fov)
            plt.show()

        assert True

    def test_GaussianDiffractionPSF_plot(self):
        psf = psfs.GaussianDiffractionPSF(diameter=32)
        fov = fovobj._centre_micado_fov()
        if PLOTS:
            psf.plot(fov)
            plt.show()

    def test_AnisocadoPsf_plot(self):
        psf = psfs.AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                                     strehl=0.5, wavelength=2.15)
        fov = fovobj._centre_micado_fov()

        if PLOTS:
            psf.plot(fov)
            plt.show()

        assert True

    def test_FieldConstantPSF_plot(self):
        psf = psfs.FieldConstantPSF(filename="test_ConstPSF.fits")
        fov = fovobj._centre_micado_fov()

        if PLOTS:
            psf.plot(fov)
            plt.show()

        assert True

# *** Not working very well for now  ***
#    def test_FieldVaryingPSF_plot(self):
#
#        psf = psfs.FieldVaryingPSF(filename="test_circular_fvpsf.fits")
#        centre_fov = _centre_fov(n=62)
#        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
#        centre_fov.hdu.data = np.zeros((nax2, nax1))
#
#        x, y = np.random.randint(6, nax1 - 6, (2, 150))
#        centre_fov.hdu.data[x, y] = 1
#    # centre_fov.hdu.data[6:nax1-6:10, 6:nax1-6:10] = 1
#        centre_fov.fields = [1]
#        psf = psf.apply_to(fov)
#        if PLOTS:
#            psf.plot(centre_fov)
#            plt.show()

#        assert True


class TestDetectorEffectsPlots:
    def test_PoorMansHxRGReadoutNoise_plot(self):
        from scopesim.detector import Detector
        level, dit, hw = 0.5, 10, 16
        from scopesim.optics.image_plane_utils import header_from_list_of_xy
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        eff = electronic.PoorMansHxRGReadoutNoise(noise_std=1,
                                                  n_channels=1, ndit=1)
        dtcr = Detector(hdr)

        if PLOTS:
            fig = plt.figure(figsize=(6, 3))
            plt.subplot(121)
            eff.plot(dtcr)
            plt.subplot(122)
            eff.plot_hist(dtcr)
            plt.show()

    def test_BasicReadoutNoise_plot(self):
        from scopesim.detector import Detector
        level, dit, hw = 0.5, 10, 16
        from scopesim.optics.image_plane_utils import header_from_list_of_xy
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        eff = electronic.BasicReadoutNoise(noise_std=1,
                                              n_channels=1, ndit=1)
        dtcr = Detector(hdr)

        if PLOTS:

            fig = plt.figure(figsize=(6, 3))
            plt.subplot(121)
            eff.plot(dtcr)
            plt.subplot(122)
            eff.plot_hist(dtcr)
            plt.show()

    def test_ShotNoise_plot(self):
        from scopesim.detector import Detector
        level, dit, hw = 0.5, 10, 16
        from scopesim.optics.image_plane_utils import header_from_list_of_xy
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        eff = electronic.ShotNoise()  # Does it do something?
        dtcr = Detector(hdr)

        if PLOTS:
            fig = plt.figure(figsize=(6, 3))
            plt.subplot(121)
            eff.plot(dtcr)
            plt.subplot(122)
            eff.plot_hist(dtcr)
            plt.show()

    def test_DarkCurrent_plot(self):
        from scopesim.detector import Detector
        level, dit, hw = 0.5, 10, 16
        from scopesim.optics.image_plane_utils import header_from_list_of_xy
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        eff = electronic.DarkCurrent(value=0.1, dit=10, ndit=20) # raise error for ints!
        dtcr = Detector(hdr)

        if PLOTS:
            eff.plot(dtcr)
            plt.show()

    def test_LinearityCurve_plot(self):
        eff = electronic.LinearityCurve(ndit=1, filename="test_linearity.dat")  # raise error for ints!

        if PLOTS:
            eff.plot()
            plt.xlim(0, 1e2)
            plt.show()


    def test_ReferencePixelBorder_plot(self):
        from scopesim.optics.image_plane import ImagePlane
        from scopesim.tests.mocks.py_objects.imagehdu_objects import _image_hdu_square

        implane = ImagePlane(_image_hdu_square().header)
        implane.hdu.data = np.ones(implane.hdu.data.shape)
        eff = electronic.ReferencePixelBorder(all=5, top=15)

        if PLOTS:
            eff.plot(implane)
            plt.show()

    def test_BinnedImage_plot(self):
        """
        Not implemented yet. What do we plot here?

        """
        from scopesim.detector import Detector
        from scopesim.optics.image_plane_utils import header_from_list_of_xy
        level, dit, hw = 0.5, 10, 16
        hdr = header_from_list_of_xy([-hw, hw], [-hw, hw], 1, "D")
        dtcr = Detector(hdr)
        eff = electronic.BinnedImage(bin_size=4)
        dtcr = eff.apply_to(dtcr)
        if PLOTS:
            plt.imshow(dtcr.data)
            plt.show()


class TestADC:

    def test_ADC_plot(self):
        atmo_params = {"airmass": 1.14,
                       "temperature": 7,
                       "humidity": 0.5,
                       "pressure": 0.755,
                       "latitude": -26,
                       "altitude": 2400,
                       "pupil_angle": 0,
                       "pixel_scale": 1,
                       "wave_min": 0.5,
                       "wave_mid": 0.5,
                       "wave_max": 2.5}

        pixel_scale = 0.04
        waves = (0.7, 0.8)
        adc = shifts.AtmosphericDispersionCorrection(**atmo_params)
        fov = _centre_fov(n=10, waverange=waves)
        fov.header["CDELT1"] = 1 / 3600 * pixel_scale
        fov.header["CDELT2"] = 1 / 3600 * pixel_scale
        old_crpix_d = np.array([fov.header["CRPIX1D"], fov.header["CRPIX2D"]])
        adc.apply_to(fov)
        new_crpix_d = np.array([fov.header["CRPIX1D"], fov.header["CRPIX2D"]])
        fov_shifts = new_crpix_d - old_crpix_d
        adc_x_shift = fov_shifts[0] * fov.header["CDELT1"] * 3600
        adc_y_shift = fov_shifts[1] * fov.header["CDELT1"] * 3600

#  Maybe create point sources and see them shift
        if PLOTS:

            print(fov_shifts)
