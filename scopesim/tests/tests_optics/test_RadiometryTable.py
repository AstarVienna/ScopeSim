# 1 read in the tables
# 2 read in curves from the set of unique files
# 3 create a dictionary of curves
#
import pytest
import os

import numpy as np
from astropy.table import Table
from astropy.io import ascii as ioascii
from astropy import units as u

from synphot import SpectralElement, SourceSpectrum

import scopesim.optics.radiometry_utils as rad_utils
from scopesim import utils
from scopesim.optics import radiometry as opt_rad
from scopesim.optics import surface as opt_surf
import scopesim as sim


MOCK_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__.insert(0, MOCK_DIR)


def synphot_version():
    from synphot.version import version
    nums = version.split(".")
    return float(nums[0]) + float(nums[1]) * 0.1


@pytest.fixture(scope="module")
def input_tables():
    filenames = ["LIST_mirrors_ELT.tbl",
                 "LIST_mirrors_SCAO_relay.tbl",
                 "LIST_mirrors_MICADO_Wide.tbl"]

    return [os.path.join(MOCK_DIR, fname) for fname in filenames]


@pytest.mark.usefixtures("input_tables")
class TestRadiometryTableInit:
    def test_initialises_with_no_input(self):
        rt = opt_rad.RadiometryTable()
        assert isinstance(rt, opt_rad.RadiometryTable)
        assert rt.table is None

    def test_initialises_with_single_table(self, input_tables):
        rt = opt_rad.RadiometryTable([input_tables[0]])
        assert isinstance(rt, opt_rad.RadiometryTable)
        assert len(rt.table) == 5

    def test_initialises_with_list_of_tables(self, input_tables):
        rt = opt_rad.RadiometryTable(input_tables)
        assert isinstance(rt, opt_rad.RadiometryTable)
        assert len(rt.table) == 19


@pytest.mark.usefixtures("input_tables")
class TestRadiometryTableAddSurfaceList:
    def test_append_single_table_from_filename(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list(input_tables[0])
        assert len(rad_table.table) == 5

    def test_combine_two_tables_from_filename(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list([input_tables[0], input_tables[1]])
        assert len(rad_table.table) == 6

    def test_combine_list_of_filename(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list(input_tables)
        assert len(rad_table.table) == 19

    def test_all_surfaces_where_added_to_dict_surface(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list(input_tables)
        names = rad_table.table["name"]
        assert np.all(name in rad_table.surface for name in names)


@pytest.mark.usefixtures("input_tables")
class TestRadiometryTableAddSurface:
    @pytest.mark.parametrize("position", (0, 2, 5))
    def test_add_empty_surface_to_full_table(self, input_tables, position):
        rt = opt_rad.RadiometryTable(tables=(input_tables[0]))
        surf = opt_surf.SpectralSurface()
        rt.add_surface(surf, "new_surf", position)
        colname = utils.real_colname("name", rt.table.colnames)
        assert rt.table[colname][position] == "new_surf"
        assert "new_surf" in rt.surfaces
        assert isinstance(rt.surfaces["new_surf"], opt_surf.SpectralSurface)

    def test_add_empty_surface_to_empty_table(self):
        rt = opt_rad.RadiometryTable()
        surf = opt_surf.SpectralSurface()
        with pytest.raises(ValueError):
            rt.add_surface(surf, "new_surf", 0)


@pytest.mark.usefixtures("input_tables")
class TestRadiometryTableGetThroughput:
    def test_return_spectral_element_from_get_throughput(self, input_tables):
        rt = opt_rad.RadiometryTable(input_tables)
        thru = rt.get_throughput()
        assert isinstance(thru, SpectralElement)

    def test_return_spectral_element_for_only_2_rows(self, input_tables):
        rt = opt_rad.RadiometryTable(input_tables)
        thru = rt.get_throughput(start=1, end=3)
        assert isinstance(thru, SpectralElement)

        if float(synphot_version()) < 0.2:
            assert thru.model.n_submodels() == 2
        else:
            assert thru.model.n_submodels == 2

    def test_return_none_for_empty_radiometry_table(self):
        rt = opt_rad.RadiometryTable()
        thru = rt.get_throughput()
        assert thru is None


@pytest.mark.usefixtures("input_tables")
class TestRadiometryTableGetEmission:
    def test_return_spectral_element_from_get_throughput(self, input_tables):
        rt = opt_rad.RadiometryTable(input_tables)
        etendue = (996*u.m**2) * (0.004*u.arcsec)**2
        emiss = rt.get_emission(etendue=etendue)
        assert isinstance(emiss, SourceSpectrum)

    def test_return_spectral_element_for_only_2_rows(self, input_tables):
        rt = opt_rad.RadiometryTable(input_tables)
        etendue = (996*u.m ** 2) * (0.004 * u.arcsec) ** 2
        emiss = rt.get_emission(etendue=etendue, start=1, end=3)
        assert isinstance(emiss, SourceSpectrum)

        if float(synphot_version()) < 0.2:
            assert emiss.model.n_submodels() == 9
        else:
            assert emiss.model.n_submodels == 7

    def test_return_none_for_empty_radiometry_table(self):
        rt = opt_rad.RadiometryTable()
        emiss = rt.get_throughput()
        assert emiss is None


@pytest.mark.usefixtures("input_tables")
class TestCombineEmissions:
    def test_super_simple_case(self):
        n = 11
        surf = opt_surf.SpectralSurface(wavelength=np.linspace(1, 2, n) * u.um,
                                        transmission=0.5*np.ones(n),
                                        area=2*u.m**2,
                                        angle=0*u.deg)
        dic = {"surf" + str(i + 1): surf for i in range(3)}
        tbl = ioascii.read(""" name action
                surf1 reflection
                surf2 transmission
                surf3 transmission """)
        etendue = 1 * u.m**2 * u.arcsec**2
        combi = rad_utils.combine_emissions(tbl, dic, [0, 1, 2], etendue)
        assert isinstance(combi, SourceSpectrum)

    def test_returns_source_spectrum_for_full_path(self, input_tables):
        rt = opt_rad.RadiometryTable(tables=input_tables)
        row_list = np.arange(len(rt.table))
        etendue = (996*u.m**2) * (0.004*u.arcsec)**2
        comb_emission = rad_utils.combine_emissions(rt.table, rt.surfaces,
                                                    row_indexes=row_list,
                                                    etendue=etendue)
        assert isinstance(comb_emission, SourceSpectrum)


@pytest.mark.usefixtures("input_tables")
class TestCombineThroughputs:
    def test_returns_spectral_element_containing_everything(self, input_tables):
        rt = opt_rad.RadiometryTable(tables=(input_tables))
        row_list = np.arange(len(rt.table))
        comb_throughput = rad_utils.combine_throughputs(rt.table, rt.surfaces,
                                                        rows_indexes=row_list)
        assert isinstance(comb_throughput, SpectralElement)

    def test_super_simple_combine_3_surfaces(self):
        n = 10
        surf = opt_surf.SpectralSurface(wavelength=np.linspace(1, 2, n)*u.um,
                                        transmission=np.ones(n))
        dic = {"surf"+str(i+1): surf for i in range(3)}
        tbl = ioascii.read(""" name action
        surf1 reflection
        surf2 transmission
        surf3 transmission """)
        combi = rad_utils.combine_throughputs(tbl, dic, [0, 1, 2])
        assert isinstance(combi, SpectralElement)

    def test_return_none_if_table_is_empty(self):
        n = 10
        surf = opt_surf.SpectralSurface(wavelength=np.linspace(1, 2, n) * u.um,
                                        transmission=np.ones(n))
        dic = {"surf" + str(i + 1): surf for i in range(3)}
        tbl = Table()
        combi = rad_utils.combine_throughputs(tbl, dic, [0, 1, 2])
        assert combi is None


@pytest.mark.usefixtures("input_tables")
class TestCombineTables:
    def test_adds_two_tables(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblB = Table(names=["colA", "colB"], data=[[2, 3], [2, 3]])
        tblC = rad_utils.combine_tables(tblB, tblA)
        assert np.all(tblC["colB"] == np.array([0, 1, 2, 3]))

    def test_adds_single_table(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblC = rad_utils.combine_tables(tblA)
        assert np.all(tblC["colA"] == np.array([0, 1]))

    def test_adds_three_tables_to_old_table(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblB = Table(names=["colA", "colB"], data=[[2, 3], [2, 3]])
        tblC = Table(names=["colA", "colB"], data=[[4, 5], [4, 5]])
        tblD = Table(names=["colA", "colB"], data=[[6, 7], [6, 7]])
        tblE = rad_utils.combine_tables([tblB, tblC, tblD], tblA)
        assert np.all(tblE["colA"] == np.arange(8))

    def test_adds_table_from_filename_to_nothing(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        tblC = rad_utils.combine_tables(tblA)
        assert len(tblC) == 5

    def test_adds_table_from_filename_to_table_object(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        tblB = input_tables[1]
        tblC = rad_utils.combine_tables(tblB, tblA)
        assert len(tblC) == 6

    def test_adds_table_from_filename_to_table_from_file(self, input_tables):
        tblA = input_tables[0]
        tblB = input_tables[1]
        tblC = rad_utils.combine_tables(tblB, tblA)
        assert len(tblC) == 6

    def test_adds_3_tables_from_filename_to_nothing(self, input_tables):
        tblC = rad_utils.combine_tables(input_tables)
        assert len(tblC) == 19

    def test_prepend_table(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblB = Table(names=["colA", "colB"], data=[[2, 3], [2, 3]])
        tblC = rad_utils.combine_tables(tblB, tblA, prepend=True)
        assert np.all(tblC["colB"] == np.array([2, 3, 0, 1]))


@pytest.mark.usefixtures("input_tables")
class TestMakeSurfaceFromRow:
    def test_return_none_from_empty_row(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        surf = rad_utils.make_surface_from_row(tblA[0])
        assert isinstance(surf, opt_surf.SpectralSurface)

    def test_surface_has_processed_ter_filename_in_row(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        surf = rad_utils.make_surface_from_row(tblA[0])
        assert isinstance(surf.transmission, SpectralElement)
        assert isinstance(surf.reflection, SpectralElement)
        assert isinstance(surf.emissivity, SpectralElement)
        assert isinstance(surf.emission, SourceSpectrum)


class TestRealColname:
    @pytest.mark.parametrize("name, colnames", [
                             ("yahoo", ["Yahoo", "Bogus"]),
                             ("yahoo", ["yahoo", "Bogus"]),
                             ("yahoo", ["YAHOO", "Bogus"])])
    def test_returns_real_name(self, name, colnames):
        assert utils.real_colname(name, colnames) == colnames[0]

    def test_returns_none_if_name_not_in_colname(self):
        assert utils.real_colname("yahoo", ["Bogus"]) is None


@pytest.mark.usefixtures("input_tables")
class TestMakeSurfaceDictFromTable:
    def test_return_dict_from_table(self, input_tables):
        tbl = ioascii.read(input_tables[0])
        surf_dict = rad_utils.make_surface_dict_from_table(tbl)
        assert isinstance(surf_dict, dict)
        assert "M1" in surf_dict


class TestInsertIntoOrderedDict:
    @pytest.mark.parametrize("dic, new_entry, pos",
                     [({},                   ["a", 1], 0),
                      ({"x": 42, "y": 3.14}, {"a": 1}, 0),
                      ({"x": 42, "y": 3.14}, ("a", 1), 2),
                      ({"x": 42, "y": 3.14}, [("b", 2), ("a", 1)], -1)])
    def test_works_as_prescribed(self, dic, new_entry, pos):
        new_dic = utils.insert_into_ordereddict(dic, new_entry, pos)
        print(new_dic, pos)
        assert list(new_dic.keys())[pos] == "a"
        assert list(new_dic.values())[pos] == 1
        assert new_dic["a"] == 1
        if "x" in new_dic:
            assert new_dic["x"] == 42
        if "b" in new_dic:
            assert new_dic["b"] == 2


class TestEmptyType:
    @pytest.mark.parametrize("x, expected",
                             [(int, 0), (float, 0.), (bool, False), (str, " ")])
    def test_works_for_all_common_types(self, x, expected):
        assert utils.empty_type(x) == expected


@pytest.mark.usefixtures("input_tables")
class TestAddSurfaceToTable:
    @pytest.mark.parametrize("position", [0, 2, 5])
    def test_(self, input_tables, position):
        tbl = ioascii.read(input_tables[0])
        surf = opt_surf.SpectralSurface(tbl[0]["filename"])
        tbl = rad_utils.add_surface_to_table(tbl, surf, "new_row", position)
        assert tbl[position]["filename"] == surf.meta["filename"]
        assert tbl[position]["name"] == "new_row"


class TestRadiometryTableFromELT:
    def local_basic_test_comparing_single_and_5_component_elt_reflections(self):
        import scopesim
        scopesim.rc.__search_path__.insert(0, ["C:\Work\irdb\ELT"])

        fname = "LIST_mirrors_ELT.tbl"
        comb = scopesim.effects.SurfaceList(filename=fname)

        fname = "TER_ELT_System_20190611.dat"
        eso = scopesim.effects.TERCurve(filename=fname)
        eso.surface.meta["temperature"] = 0

        from matplotlib import pyplot as plt
        wave = np.arange(0.3, 3, 0.001) * 1e4
        plt.plot(wave * 1E-4, eso.surface.reflection(wave), label="ESO 2019")
        plt.plot(wave * 1E-4, comb.throughput(wave), label="5 component")

        plt.xlabel("Wavelength [um]")
        plt.ylabel("Throughput ")
        plt.legend(loc=5)
        plt.xlim(0.3, 3)

        plt.show()

        print(comb.throughput(3E4), eso.surface.reflection(3E4))
