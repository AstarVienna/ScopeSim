""" Tests for the the new SurfaceList object """
from pytest import approx
import os
import numpy as np
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt

from synphot import SourceSpectrum, SpectralElement

from scopesim.effects import surface_list as sl
from scopesim.optics import SpectralSurface

from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim import rc


MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]

PLOTS = False


def micado_surf_list():
    return sl.SurfaceList(filename="LIST_mirrors_MICADO_Wide.tbl")


def surf_list_kwargs(n=11):
    kwargs = {"array_dict": {"name": ["M{}".format(i) for i in range(n)],
                             "area": [1.0] * n,
                             "angle": [0] * n,
                             "temperature": [0] * n,
                             "action": ["reflection"] * n,
                             "filename": ["TER_mirror_gold.dat"] * n},
              "outer_unit": "m",
              "inner_unit": "m",
              "angle_unit": "deg",
              "temperature_unit": "deg_C",
              "etendue": (1 * u.m * u.arcsec) ** 2}
    return kwargs


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(sl.SurfaceList(), sl.SurfaceList)

    def test_initialises_with_valid_filename(self):
        surf_list = micado_surf_list()
        assert isinstance(surf_list, sl.SurfaceList)
        assert isinstance(surf_list.data, Table)
        for key in surf_list.surfaces:
            assert isinstance(surf_list.surfaces[key], SpectralSurface)

    def test_initialises_from_array_dict(self):
        n = 5
        kwargs = surf_list_kwargs(n)
        surf_list = sl.SurfaceList(**kwargs)
        assert len(surf_list.data) == n
        assert isinstance(surf_list, sl.SurfaceList)
        assert surf_list.surfaces["M1"].area == 1 * u.m**2


class TestGetEmission:
    def test_returns_source_spectrum_object(self):
        surf_list = micado_surf_list()
        etendue = (996 * u.m ** 2) * (0.004 * u.arcsec) ** 2
        assert isinstance(surf_list.get_emission(etendue), SourceSpectrum)

        if PLOTS:
            wave = np.arange(10, 200, 1) * u.um
            plt.semilogy(wave, surf_list.get_emission(etendue)(wave))
            plt.show()

    def test_combines_emission_correctly(self):
        n = 5
        kwargs = surf_list_kwargs(n)
        surf_list = sl.SurfaceList(**kwargs)

        wave = np.arange(0.8, 2.5, 0.01) * u.um
        sum_values = np.sum([surf_list.surfaces[key].emission(wave).value
                             for key in surf_list.surfaces], axis=0)
        comb_value = surf_list.emission(wave).value
        sollwert = np.sum([0.985 ** n for n in range(n)]) / n

        assert np.average(comb_value / sum_values) == approx(sollwert, rel=1e-5)

        # print(sum_values, comb_value, comb_value / sum_values, sollwert)
        if PLOTS:
            wave = np.linspace(0.8, 2.5, 100) * u.um
            for key, surf in surf_list.surfaces.items():
                plt.semilogy(wave, surf.emission(wave))
            plt.semilogy(wave, surf_list.emission(wave))
            plt.show()


class TestGetThroughput:
    def test_combines_throughputs(self):
        surf_list = micado_surf_list()
        av_value = surf_list.throughput([2] * u.um)[0].value
        assert av_value == approx(0.985**11, rel=0.01)
        assert isinstance(surf_list.throughput, SpectralElement)

        if PLOTS:
            wave = np.linspace(0.8, 2.5) * u.um
            for key, surf in surf_list.surfaces.items():
                plt.plot(wave, surf.throughput(wave))
            plt.plot(wave, surf_list.throughput(wave))
            plt.show()


class TestApplyTo:
    def test_adds_bg_to_source_if_source_has_no_bg(self):
        n = 11
        kwargs = surf_list_kwargs(n)
        surf_list = sl.SurfaceList(**kwargs)
        src = so._vega_source()

        wave = np.linspace(1, 2.5, 1000) * u.um
        flux_before = src.spectra[0](wave).value
        src = surf_list.apply_to(src)
        flux_after = src.spectra[0](wave).value

        assert np.average(flux_after / flux_before) == approx(0.985**n)
        assert src.fields[-1].header["BG_SRC"] is True
        assert src.fields[-1].header["SPEC_REF"] == 1

        if PLOTS:
            plt.semilogy(wave, flux_before)
            plt.semilogy(wave, flux_after)
            plt.show()


class TestPlot:
    def test_plotting(self):
        if PLOTS:
            surf_list = micado_surf_list()
            surf_list = sl.SurfaceList(**surf_list_kwargs(2))
            surf_list.plot("xe")
            plt.show()
