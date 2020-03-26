import pytest
from pytest import approx

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

import scopesim.optics.fov_utils
from scopesim.optics import fov

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scopesim.tests.mocks.py_objects.source_objects import _table_source, \
    _image_source, _combined_source
from scopesim.tests.mocks.py_objects.header_objects import _basic_fov_header


PLOTS = False


@pytest.fixture(scope="function")
def table_source():
    return _table_source()


@pytest.fixture(scope="function")
def image_source():
    return _image_source()


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
            fov.FieldOfView()

    def test_initialises_with_header_and_waverange(self, basic_fov_header):
        print(dict(basic_fov_header))
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        assert isinstance(the_fov, fov.FieldOfView)

    def test_throws_error_if_no_wcs_in_header(self):
        with pytest.raises(ValueError):
            fov.FieldOfView(fits.Header(), (1, 2) * u.um, area=1*u.m**2)


@pytest.mark.usefixtures("basic_fov_header", "table_source", "image_source")
class TestFieldOfViewExtractFrom:
    def test_creates_single_combined_from_multiple_tables(self, table_source,
                                                          basic_fov_header):
        src = table_source + table_source
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 1
        assert isinstance(the_fov.fields[0], Table)

    def test_creates_single_combined_from_multiple_images(self, image_source,
                                                          basic_fov_header):
        src = image_source + image_source
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 1
        assert isinstance(the_fov.fields[0], fits.ImageHDU)

    def test_creates_two_fields_for_tables_and_images(self, basic_fov_header,
                                                      image_source,
                                                      table_source):
        src = image_source + table_source
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 2
        assert isinstance(the_fov.fields[0], Table)
        assert isinstance(the_fov.fields[1], fits.ImageHDU)

    def test_ignores_fields_outside_fov_boundary(self, basic_fov_header):

        src = _combined_source(dx=[200, 200, 200])
        src.fields[0]["x"] += 200

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)

        assert len(the_fov.fields) == 0


@pytest.mark.usefixtures("basic_fov_header")
class TestFieldOfViewView:
    def test_views_with_only_image(self, basic_fov_header):
        src = _image_source()
        flux = src.photons_in_range(1*u.um, 2*u.um).value
        orig_sum = np.sum(src.fields[0].data * flux)

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
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
        src = _image_source(angle=30)
        flux = src.photons_in_range(1*u.um, 2*u.um).value
        orig_sum = np.sum(src.fields[0].data) * flux

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.isclose(np.sum(view), orig_sum, rtol=1e-3)

        if PLOTS:
            plt.imshow(src.fields[0].data, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

            plt.imshow(view, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_only_table(self, basic_fov_header):
        src = _table_source()
        fluxes = src.photons_in_range(1*u.um, 2*u.um)
        phs = fluxes[src.fields[0]["ref"]] * src.fields[0]["weight"]

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2) * u.um, area=1*u.m**2)
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

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um, area=1*u.m**2)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.sum(the_fov.fields[0]["flux"]) == approx(36)
        assert np.isclose(np.sum(the_fov.fields[1].data), orig_sum)
        assert np.isclose(np.sum(view), orig_sum + 36)

        if PLOTS:
            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()


@pytest.mark.usefixtures("basic_fov_header", "image_source", "table_source")
class TestIsFieldInFOV:
    def test_returns_true_for_table_inside(self, basic_fov_header,
                                           table_source):
        assert scopesim.optics.fov_utils.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_false_for_table_outside(self, basic_fov_header,
                                             table_source):
        table_source.fields[0]["x"] += 200
        assert not scopesim.optics.fov_utils.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_true_for_table_corner_inside(self, basic_fov_header,
                                                  table_source):
        table_source.fields[0]["x"] += 9
        table_source.fields[0]["y"] += 14
        assert scopesim.optics.fov_utils.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_true_for_image_inside(self, basic_fov_header,
                                           image_source):
        assert scopesim.optics.fov_utils.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_false_for_image_outside(self, basic_fov_header,
                                             image_source):
        image_source.fields[0].header["CRVAL1"] += 200*u.arcsec.to(u.deg)
        assert not scopesim.optics.fov_utils.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_true_for_image_in_corner(self, basic_fov_header,
                                              image_source):
        image_source.fields[0].header["CRVAL1"] += 10*u.arcsec.to(u.deg)
        image_source.fields[0].header["CRVAL2"] -= 10*u.arcsec.to(u.deg)
        assert scopesim.optics.fov_utils.is_field_in_fov(basic_fov_header, image_source.fields[0])


class TestMakeFluxTable:
    def test_flux_in_equals_flux_out(self):
        pass


class TestCombineTableFields:
    def test_flux_in_equals_flux_out(self):
        pass


class TestCombineImageHDUFields:
    def test_flux_in_equals_flux_out(self):
        pass
