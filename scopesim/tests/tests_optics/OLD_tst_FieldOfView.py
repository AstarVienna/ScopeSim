import pytest
from pytest import approx

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from scopesim.optics.fov import FieldOfView
from scopesim.optics import fov_utils

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scopesim.tests.mocks.py_objects.source_objects import _table_source, \
    _image_source, _combined_source, _cube_source
from scopesim.tests.mocks.py_objects.header_objects import _basic_fov_header


PLOTS = False


@pytest.fixture(scope="function")
def table_source():
    return _table_source()


@pytest.fixture(scope="function")
def image_source():
    return _image_source()


@pytest.fixture(scope="function")
def cube_source():
    return _cube_source()


@pytest.fixture(scope="function")
def combined_source():
    return _combined_source()


@pytest.fixture(scope="function")
def basic_fov_header():
    return _basic_fov_header()


@pytest.mark.usefixtures("basic_fov_header")
class TestFieldOfViewInit:
    def test_initialises_with_nothing_raise_error(self):
        with pytest.raises(TypeError):
            FieldOfView()

    def test_initialises_with_header_and_waverange(self, basic_fov_header):
        print(dict(basic_fov_header))
        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        assert isinstance(the_fov, FieldOfView)

    def test_throws_error_if_no_wcs_in_header(self):
        with pytest.raises(ValueError):
            FieldOfView(fits.Header(), (1, 2) * u.um, area=1*u.m**2)



@pytest.mark.usefixtures("basic_fov_header", "table_source", "image_source",
                         "cube_source")
class TestFieldOfViewExtractFrom:
    # Note: Deleted tests for combining fields - why combine Source fields if you don't need to
    def test_contains_all_fields_inside_fov(self, basic_fov_header, cube_source,
                                            image_source, table_source):
        src = image_source + cube_source + table_source
        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 3
        assert isinstance(the_fov.fields[0], fits.ImageHDU)
        assert isinstance(the_fov.fields[1], fits.ImageHDU)
        assert  the_fov.fields[1].header["NAXIS"] == 3
        assert isinstance(the_fov.fields[2], Table)

    # OLD FOV - Can't work with new FOV because new FOV has copies of spectra

    # def test_all_spectra_are_referenced_correctly(self, basic_fov_header,
    #                                               image_source, table_source,
    #                                               cube_source):
    #     src = image_source + cube_source + table_source
    #     the_fov = FieldOfView(basic_fov_header, (1, 2) * u.um,
    #                               area=1 * u.m ** 2)
    #     the_fov.extract_from(src)
    #     # check the same spectrum object is referenced by both lists
    #     assert all([the_fov.spectra[i] == src.spectra[i]
    #                 for i in the_fov.spectra])

    def test_ignores_fields_outside_fov_boundary(self, basic_fov_header):

        src = _combined_source(dx=[200, 200, 200])
        src.fields[0]["x"] += 200

        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)

        assert len(the_fov.fields) == 0

    def test_contains_cube_when_passed_a_cube_source(self, basic_fov_header,
                                                     cube_source):
        src = cube_source
        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)

        assert len(the_fov.fields) == 1
        assert len(the_fov.spectra) == 0
        assert len(the_fov.fields[0].data.shape) == 3


@pytest.mark.usefixtures("basic_fov_header")
class TestFieldOfViewView:
    def test_views_with_only_image(self, basic_fov_header):
        src = _image_source()
        flux = src.photons_in_range(1*u.um, 2*u.um).value
        orig_sum = np.sum(src.fields[0].data * flux)

        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        the_fov.view()

        assert np.isclose(np.sum(the_fov.hdu.data), orig_sum)

        if PLOTS:
            plt.imshow(src.fields[0].data, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

            plt.imshow(the_fov.hdu.data, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_rotated_image(self, basic_fov_header):
        area = 1 * u.m ** 2

        # works for angle=0, not working for angle!=0
        src = _image_source(angle=30)
        flux = src.photons_in_range(1*u.um, 2*u.um, area=area).value
        in_sum = np.sum(src.fields[0].data) * flux

        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=area)
        the_fov.extract_from(src)
        view = the_fov.view()
        out_sum = np.sum(view)

        assert out_sum == approx(in_sum, rel=1e-3)

        if PLOTS:
            plt.subplot(121)
            plt.imshow(src.fields[0].data, origin="lower", norm=LogNorm())
            plt.colorbar()

            plt.subplot(122)
            plt.imshow(view, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_only_table(self, basic_fov_header):
        src = _table_source()
        fluxes = src.photons_in_range(1*u.um, 2*u.um)
        phs = fluxes[src.fields[0]["ref"]] * src.fields[0]["weight"]

        the_fov = FieldOfView(basic_fov_header, (1, 2) * u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.sum(view) == np.sum(phs).value - 4

        if PLOTS:
            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_tables_and_images(self, basic_fov_header):
        src = _combined_source(im_angle=0, weight=[1, 1, 1],
                               dx=[-2, 0, 0], dy=[2, 2, 0])
        ii = src.fields[3].header["SPEC_REF"]
        flux = src.photons_in_range(1*u.um, 2*u.um, indexes=[ii]).value
        orig_sum = np.sum(src.fields[3].data) * flux

        the_fov = FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.isclose(np.sum(the_fov.fields[1].data), orig_sum)
        assert np.isclose(np.sum(view), orig_sum + 36)
        assert np.sum(the_fov.fields[0]["flux"]) == approx(36)

        if PLOTS:
            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()


@pytest.mark.usefixtures("basic_fov_header", "image_source", "table_source")
class TestIsFieldInFOV:
    def test_returns_true_for_table_inside(self, basic_fov_header,
                                           table_source):
        assert fov_utils.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_false_for_table_outside(self, basic_fov_header,
                                             table_source):
        table_source.fields[0]["x"] += 200
        assert not fov_utils.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_true_for_table_corner_inside(self, basic_fov_header,
                                                  table_source):
        table_source.fields[0]["x"] += 9
        table_source.fields[0]["y"] += 14
        assert fov_utils.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_true_for_image_inside(self, basic_fov_header,
                                           image_source):
        assert fov_utils.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_false_for_image_outside(self, basic_fov_header,
                                             image_source):
        image_source.fields[0].header["CRVAL1"] += 200*u.arcsec.to(u.deg)
        assert not fov_utils.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_true_for_image_in_corner(self, basic_fov_header,
                                              image_source):
        image_source.fields[0].header["CRVAL1"] += 10*u.arcsec.to(u.deg)
        image_source.fields[0].header["CRVAL2"] -= 10*u.arcsec.to(u.deg)
        assert fov_utils.is_field_in_fov(basic_fov_header, image_source.fields[0])



@pytest.mark.usefixtures("cube_source")
class TestCubeSourceInput:
    def test_source_cube_exists(self, cube_source):
        assert len(cube_source.fields[0].data.shape) == 3
