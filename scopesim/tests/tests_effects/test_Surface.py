# -*- coding: utf-8 -*-
import pytest
from pytest import approx

import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt
from synphot import SpectralElement, SourceSpectrum

from scopesim.effects import TERCurve
from scopesim.optics.surface import SpectralSurface, scale_spectrum
from scopesim.source.source_templates import vega_spectrum


PLOTS = False


class TestTERCurveInit:
    def test_initialise_with_nothing(self):
        assert isinstance(TERCurve(), TERCurve)

    def test_initalises_with_list_of_surfaces(self, mock_path_micado):
        filename = str(mock_path_micado / "TC_filter_Ks.dat")
        surf = TERCurve(filename=filename)
        assert isinstance(surf, TERCurve)

    def test_initalises_with_two_arrays(self):
        surf = TERCurve(wavelength=np.array([0.5, 2.5]),
                        transmission=np.array([1, 1]),
                        wavelength_unit="um")
        assert isinstance(surf, TERCurve)


class TestSurfaceAttribute:
    def test_returns_surface_object(self, mock_path_micado):
        filename = str(mock_path_micado / "TC_filter_Ks.dat")
        surf = TERCurve(filename=filename)

        assert isinstance(surf.surface, SpectralSurface)
        assert isinstance(surf.surface.transmission, SpectralElement)
        assert isinstance(surf.surface.emission, SourceSpectrum)

    def test_returns_surface_object_for_arrays(self):
        surf = TERCurve(wavelength=[0.5, 1.5, 2.5],
                        transmission=[0.1, 0.1, 0.1],
                        emission=[1, 2, 3]*u.ph/u.s/u.um/u.m**2,
                        wavelength_unit="um")

        assert isinstance(surf.surface, SpectralSurface)
        assert isinstance(surf.surface.transmission, SpectralElement)
        assert isinstance(surf.surface.emission, SourceSpectrum)

        if PLOTS:
            wave = surf.surface.wavelength
            plt.plot(wave, surf.surface.emission(wave))
            plt.show()


@pytest.mark.webtest
@pytest.mark.usefixtures("no_file_error")
class TestScaleSpectrum:
    def test_scales_vega_spectrum_to_vega_ab_or_jansky(self):
        spec = vega_spectrum()
        vega_185 = scale_spectrum(spec, "Ks", -1.85 * u.mag)
        ab_0 = scale_spectrum(spec, "Ks", 0 * u.ABmag)
        jy_3630 = scale_spectrum(spec, "Ks", 3630 * u.Jy)

        wave = np.linspace(1.8, 2.5, 1000) * u.um
        assert vega_185(wave).value == approx(ab_0(wave).value, rel=1e-2)
        assert vega_185(wave).value == approx(jy_3630(wave).value, rel=1e-2)

        if PLOTS:
            plt.plot(wave, spec(wave), "b")
            plt.plot(wave, vega_185(wave), "r:")
            plt.plot(wave, ab_0(wave), "g--")
            plt.plot(wave, jy_3630(wave), "y-.")
            plt.semilogy()
            plt.show()
