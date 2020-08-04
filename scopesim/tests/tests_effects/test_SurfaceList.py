""" Tests for the the new SurfaceList object """
from pytest import approx
import os
import numpy as np
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt

from synphot import SourceSpectrum, SpectralElement

import scopesim.effects
from scopesim.effects import surface_list as sl
from scopesim.optics import SpectralSurface
from scopesim import rc

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]

PLOTS = True


def micado_surf_list():
    return sl.SurfaceList(filename="LIST_mirrors_MICADO_Wide.tbl")


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
        kwargs = {"array_dict": {"name": ["I01_Fold1"],
                                 "outer": [0.5],
                                 "inner": [0.0],
                                 "angle": [45],
                                 "temperature": [-190],
                                 "action": ["reflection"],
                                 "filename": ["TER_mirror_gold.dat"]},
                  "outer_unit": "m",
                  "inner_unit": "m",
                  "angle_unit": "deg",
                  "temperature_unit": "deg_C"}
        surf_list = sl.SurfaceList(**kwargs)
        assert len(surf_list.data) == 1
        assert isinstance(surf_list, sl.SurfaceList)
        assert surf_list.surfaces["I01_Fold1"].area == (0.25 * u.m)**2 * np.pi


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
        n = 2
        kwargs = {"array_dict": {"name": ["M{}".format(i) for i in range(n)],
                                 "outer": [1.0]*n,
                                 "inner": [0.0]*n,
                                 "angle": [0]*n,
                                 "temperature": [0]*n,
                                 "action": ["reflection"]*2,
                                 "filename": ["TER_mirror_gold.dat"]*2},
                  "outer_unit": "m",
                  "inner_unit": "m",
                  "angle_unit": "deg",
                  "temperature_unit": "deg_C",
                  "etendue": (1*u.m*u.arcsec)**2}
        surf_list = sl.SurfaceList(**kwargs)
        m1_value = surf_list.surfaces["M0"].emission(2 * u.um).value
        m2_value = surf_list.surfaces["M1"].emission(2 * u.um).value
        comb_value = surf_list.emission(2 * u.um).value
        print(m1_value + m2_value, comb_value)
        raise ValueError("m1+m2 < comb!!!")

        if not PLOTS:
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

