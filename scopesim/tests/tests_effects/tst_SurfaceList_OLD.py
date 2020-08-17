import pytest
import os
from copy import deepcopy

import numpy as np
from astropy import units as u

from synphot import SourceSpectrum

from scopesim import rc
from scopesim.effects.surface_list_OLD import SurfaceList
from scopesim.optics.radiometry import RadiometryTable

from scopesim.tests.mocks.py_objects.source_objects import _image_source
from scopesim.tests.mocks.py_objects.effects_objects import _surf_list, \
    _surf_list_empty, _filter_surface

from matplotlib import pyplot as plt

PLOTS = False

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]



@pytest.fixture(scope="function")
def surf_list():
    return _surf_list()


@pytest.fixture(scope="function")
def surf_list_empty():
    return _surf_list_empty()


@pytest.fixture(scope="function")
def filter_surface():
    return _filter_surface()


@pytest.fixture(scope="function")
def image_source():
    return _image_source()


@pytest.mark.usefixtures("surf_list")
class TestInit:
    def test_initialise_with_nothing(self):
        assert isinstance(SurfaceList(), SurfaceList)

    def test_initalises_with_list_of_surfaces(self, surf_list):
        assert isinstance(surf_list, SurfaceList)

    def test_initialises_with_array_dict_of_surfaces(self):
        surf_list = SurfaceList(array_dict={"name": ["M1"],
                                            "outer": [1],
                                            "angle": [0],
                                            "temperature": [0],
                                            "action": ["transmission"],
                                            "wavelength": [[0.8, 2.5]],
                                            "transmission": [[0.1, 1]]
                                            },
                                outer_unit="m",
                                angle_unit="deg",
                                temperature_unit="deg_c",
                                wavelength_unit="um")
        wave = np.arange(0.8, 2.5, 0.01) * u.um
        assert np.all(surf_list.throughput(wave) > 0)
        assert isinstance(surf_list, SurfaceList)


@pytest.mark.usefixtures("surf_list")
class TestRadiometryTableAttribute:
    def test_returns_radiometry_table_object(self, surf_list):
        assert isinstance(surf_list.radiometry_table, RadiometryTable)
        assert isinstance(surf_list.get_emission(), SourceSpectrum)

        if PLOTS:
            wave = np.arange(10, 200, 1) * u.um
            plt.plot(wave, surf_list.get_emission()(wave))
            plt.semilogy()
            plt.show()


@pytest.mark.usefixtures("surf_list")
class TestAddSurfaceList:
    def test_second_list_is_joined_to_first(self, surf_list):
        surf_list_copy = deepcopy(surf_list)
        len1 = len(surf_list.radiometry_table.table)
        len2 = len(surf_list_copy.radiometry_table.table)

        surf_list.add_surface_list(surf_list_copy, prepend=True)
        len3 = len(surf_list.radiometry_table.table)

        assert len3 == len2 + len1


@pytest.mark.usefixtures("surf_list", "filter_surface")
class TestAddSurface:
    def test_extra_surface_is_joined_to_list(self, surf_list, filter_surface):
        len1 = len(surf_list.radiometry_table.table)
        surf_list.add_surface(filter_surface, "filter")

        assert len(surf_list.radiometry_table.table) == len1 + 1

        if PLOTS:
            wave = np.linspace(0.5, 2.5, 100)*u.um
            plt.subplot(121)
            plt.plot(wave, surf_list.get_throughput()(wave))

            plt.subplot(122)
            plt.plot(wave, surf_list.get_emission()(wave))
            plt.semilogy()
            plt.show()


@pytest.mark.usefixtures("surf_list", "filter_surface")
class TestFovGrid:
    def test_finds_borders_of_filter(self, filter_surface, surf_list):
        surf_list.add_surface(filter_surface, "filter")
        surf_list.meta["minimum_throughput"] = 2e-4
        waverange = surf_list.fov_grid("waveset",
                                       wave_min=0.5*u.um, wave_max=2.5*u.um)

        # assuming surf is the K-filter
        assert waverange[0] > 1.9*u.um
        assert waverange[1] < 2.4*u.um


@pytest.mark.usefixtures("surf_list_empty", "filter_surface", "image_source")
class TestApplyTo:
    def test_applied_to_source(self, filter_surface, surf_list_empty, image_source):
        surf_list_empty.add_surface(filter_surface, "filter")
        new_source = surf_list_empty.apply_to(deepcopy(image_source))

        wave = np.linspace(1.5, 2.5, 100)*u.um
        old_flux = image_source.spectra[0](wave)
        new_flux = new_source.spectra[0](wave)

        assert np.all(new_flux < old_flux)

        if PLOTS:
            plt.plot(wave, image_source.spectra[0](wave))
            plt.plot(wave, new_source.spectra[0](wave))
            plt.show()
