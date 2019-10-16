# datacontainer must read in the data if it is an ASCII file, or open a file
# handle to it if it is a FITS file
# the header(s) must be accessible as dictionaries
# if the data is in table format, the table command accesses this
# if the data is in image format, the image command accesses this

import os
import inspect
import pytest

import numpy as np
from astropy.table import Table
from astropy import units as u

from scopesim.effects.data_container import DataContainer


def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "../mocks/MICADO_SCAO_WIDE/"

    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))


MOCK_DIR = mock_dir()


@pytest.fixture(scope="module")
def data_files():
    filenames = ["PSF_basic.fits", "TC_filter_Ks.dat"]
    abs_paths = [os.path.join(MOCK_DIR, fname) for fname in filenames]

    return abs_paths


@pytest.mark.usefixtures("data_files")
class TestInit:
    def test_initialised_with_no_input(self):
        dat = DataContainer()
        assert isinstance(dat, DataContainer)

    def test_initialised_with_psf_input(self, data_files):
        dat = DataContainer(data_files[0])
        assert isinstance(dat, DataContainer)
        assert dat.is_fits is True

    def test_initialised_with_ascii_input(self, data_files):
        dat = DataContainer(data_files[1])
        assert isinstance(dat, DataContainer)
        assert dat.is_fits is False

    def test_initialised_with_arrays_dict_input(self):
        array_dict = {"wavelength" : np.linspace(1, 2, 11)*u.um,
                      "transmission" : np.ones(11)}
        dat = DataContainer(array_dict=array_dict)
        assert isinstance(dat, DataContainer)
        assert dat.is_fits is False


@pytest.mark.usefixtures("data_files")
class TestGetData:
    def test_no_file_returns_no_input(self):
        dat = DataContainer()
        assert dat.get_data() is None

    @pytest.mark.parametrize("ext, layer, datatype, dims",
                             [(0, None, type(None), None),
                              (1, None, Table, None),
                              (2, None, np.ndarray, 3),
                              (2, 1, np.ndarray, 2)])
    def test_psf_input_returns_arrays_or_tables(self, data_files, ext,
                                                layer, datatype, dims):
        datc = DataContainer(data_files[0])
        data = datc.get_data(ext, layer)
        if dims is not None:
            assert len(data.shape) == dims
        assert isinstance(data, datatype)

    def test_ascii_input_returns_table(self, data_files):
        datc = DataContainer(data_files[1])
        data = datc.get_data()
        assert isinstance(data, Table)

    def test_array_input_returns_table(self):
        array_dict = {"wavelength" : np.linspace(1, 2, 11)*u.um,
                      "transmission" : np.ones(11)}
        datc = DataContainer(array_dict=array_dict)
        data = datc.get_data()
        assert isinstance(data, Table)

