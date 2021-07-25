import pytest
from matplotlib import pyplot as plt
from scopesim.optics import FieldOfView, fov_utils
from scopesim.optics import image_plane_utils as imp_utils

from scopesim.tests.mocks.py_objects.header_objects import _basic_fov_header
from scopesim.tests.mocks.py_objects.source_objects import _cube_source

PLOTS = True


@pytest.fixture(scope="function")
def cube_source():
    return _cube_source()


@pytest.fixture(scope="function")
def basic_fov_header():
    return _basic_fov_header()


@pytest.mark.usefixtures("cube_source", "basic_fov_header")
class TestExtractAreaFromImageHDU:
    def test_returns_full_cube_for_thick_fov(self, cube_source,
                                             basic_fov_header):
        fov = FieldOfView(basic_fov_header, [0.3, 2.7])
        field = cube_source.fields[0]
        new_field = fov_utils.extract_area_from_imagehdu(field, fov.volume())

        if PLOTS:
            x, y = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(x, y, c="r")
            x, y = imp_utils.calc_footprint(field.header)
            plt.fill(x, y, c="y")
            x, y = imp_utils.calc_footprint(new_field.header)
            plt.fill(x, y, c="g")

            plt.show()

        assert new_field.header["NAXIS1"] == field.header["NAXIS1"]
        assert new_field.header["NAXIS2"] == field.header["NAXIS2"]
        assert new_field.header["NAXIS3"] == field.header["NAXIS3"]

    def test_returns_wavelength_cut_cube_for_thin_fov(self, cube_source,
                                                      basic_fov_header):
        fov = FieldOfView(basic_fov_header, [1.3, 1.7])
        field = cube_source.fields[0]
        new_field = fov_utils.extract_area_from_imagehdu(field, fov.volume())

        if PLOTS:
            x, y = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(x, y, c="r")
            x, y = imp_utils.calc_footprint(field.header)
            plt.fill(x, y, c="y")
            x, y = imp_utils.calc_footprint(new_field.header)
            plt.fill(x, y, c="g")

            plt.show()

        assert new_field.header["NAXIS1"] == field.header["NAXIS1"]
        assert new_field.header["NAXIS2"] == field.header["NAXIS2"]
        assert new_field.header["NAXIS3"] == 20

    def test_returns_eigth_cube_for_3d_offset_fov(self, cube_source,
                                                         basic_fov_header):
        hdr = basic_fov_header
        hdr["CRVAL1"] += 75 * hdr["CDELT1"]
        hdr["CRVAL2"] += 75 * hdr["CDELT2"]
        fov = FieldOfView(hdr, [1.5, 2.7])
        field = cube_source.fields[0]
        new_field = fov_utils.extract_area_from_imagehdu(field, fov.volume())

        if PLOTS:
            x, y = imp_utils.calc_footprint(basic_fov_header)
            plt.fill(x, y, c="r")
            x, y = imp_utils.calc_footprint(field.header)
            plt.fill(x, y, c="y")
            x, y = imp_utils.calc_footprint(new_field.header)
            plt.fill(x, y, c="g")

            plt.show()

        assert new_field.header["NAXIS1"] == 26
        assert new_field.header["NAXIS2"] == 26
        assert new_field.header["NAXIS3"] == 51
