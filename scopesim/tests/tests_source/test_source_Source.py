# -*- coding: utf-8 -*-

import pytest
from pytest import approx

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table
from astropy import units as u
from astropy import wcs

from synphot import SourceSpectrum
from synphot.models import Empirical1D
from synphot.units import PHOTLAM

from scopesim.source import source_utils
from scopesim.source.source import Source
from scopesim.source.source_fields import CubeSourceField, TableSourceField

from scopesim.utils import convert_table_comments_to_dict

from scopesim.tests.mocks.py_objects import source_objects as so


PLOTS = False


@pytest.fixture(scope="module")
def input_files(mock_path):
    filenames = ["test_image.fits", "test_table.fits", "test_table.tbl",
                 "test_spectrum_Flam.dat", "test_spectrum_photlam.dat"]
    filenames = [str(mock_path / fname) for fname in filenames]
    return filenames


@pytest.fixture(scope="module")
def input_hdulist(mock_path):
    filenames = ["test_image.fits"]
    filenames = [str(mock_path / fname) for fname in filenames]
    hdu_handle = fits.open(filenames[0])

    return hdu_handle


@pytest.fixture(scope="module")
def input_tables(mock_path):
    filenames = ["test_table.fits", "test_table.tbl"]
    filenames = [str(mock_path / fname) for fname in filenames]
    tbls = []
    tbls += [Table.read(filenames[0])]
    tbls += [Table.read(filenames[1], format="ascii.basic")]
    tbls[1].meta.update(convert_table_comments_to_dict(tbls[1]))

    return tbls


@pytest.fixture(scope="module")
def input_spectra(mock_path):
    filenames = ["test_spectrum_photlam.dat", "test_spectrum_Flam.dat"]
    filenames = [str(mock_path / fname) for fname in filenames]
    tbls = [ioascii.read(fname) for fname in filenames]
    specs = []
    for tbl in tbls:
        tbl.meta.update(convert_table_comments_to_dict(tbl))
        wave = tbl["wavelength"] * u.Unit(tbl.meta["wavelength_unit"])
        flux = tbl["flux"] * u.Unit(tbl.meta["flux_unit"])
        specs += [SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)]

    return specs


@pytest.fixture(scope="function")
def table_source():
    n = 100
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=4 * np.ones(n) * unit),
             SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit),
             SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n)[::-1] * unit)]
    tbl = Table(names=["x", "y", "ref", "weight"],
                data=[[5,  0, -5,  0]*u.arcsec,
                      [5, -10, 5,  0] * u.arcsec,
                      [2,  0,  1,  0],
                      [1,  1,  1,  2]])
    tbl_source = Source(table=tbl, spectra=specs)

    return tbl_source


@pytest.fixture(scope="function")
def image_source():
    n = 50
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit)]

    n = 51
    im_wcs = wcs.WCS(naxis=2)
    im_wcs.wcs.cunit = [u.arcsec, u.arcsec]
    im_wcs.wcs.cdelt = [0.2, 0.2]
    im_wcs.wcs.crval = [0, 0]
    im_wcs.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    # im_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    im_wcs.wcs.ctype = ["LINEAR", "LINEAR"]

    im = np.ones((n, n)) * 1E-11
    im[0, n-1] += 5
    im[n-1, 0] += 5
    im[n//2, n//2] += 10

    im_hdu = fits.ImageHDU(data=im, header=im_wcs.to_header())
    im_hdu.header["SPEC_REF"] = 0
    im_source = Source(image_hdu=im_hdu, spectra=specs)

    return im_source


class TestSourceInit:
    def test_initialises_with_nothing(self):
        src = Source()
        assert isinstance(src, Source)
        src.shift(0.1, 0.2)

    @pytest.mark.parametrize("ii", [0, 1])
    def test_initialises_with_table_and_2_spectrum(self, ii,
                                                   input_tables,
                                                   input_spectra):
        table = input_tables[ii]
        src = Source(table=table, spectra=input_spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, Table)
        src.shift(0.1, 0.2)

    def test_initialises_with_image_and_1_spectrum(self, input_hdulist,
                                                   input_spectra):
        src = Source(image_hdu=input_hdulist[0], spectra=input_spectra[0])
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, fits.ImageHDU)
        src.shift(0.1, 0.2)

    def test_initialises_with_image_and_flux(self, input_hdulist):
        src = Source(image_hdu=input_hdulist[0], flux=20*u.ABmag)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, fits.ImageHDU)
        src.shift(0.1, 0.2)

    def test_initialises_with_only_image(self, input_hdulist):
        input_hdulist[0].header["BUNIT"] = "ph s-1 cm-2 AA-1"
        src = Source(image_hdu=input_hdulist[0])
        assert len(src.spectra) == 1
        assert src.fields[0].header["SPEC_REF"] == 0
        src.shift(0.1, 0.2)

    def test_initialises_with_only_imagehdu_and_arcsec2(self):
        hdu = fits.ImageHDU(data=np.ones([3, 3]))
        hdu.header["BUNIT"] = "Jy/arcsec2"
        hdu.header["CRVAL1"] = 0.0
        hdu.header["CDELT1"] = 0.1
        hdu.header["CUNIT1"] = "arcsec"
        hdu.header["CRVAL2"] = 0.0
        hdu.header["CDELT2"] = 0.1
        hdu.header["CUNIT2"] = "arcsec"
        src = Source(image_hdu=hdu)

        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, fits.ImageHDU)
        src.shift(0.1, 0.2)

    @pytest.mark.parametrize("ii, dtype",
                             [(0, fits.ImageHDU),
                              (1, Table),
                              (2, Table)])
    def test_initialises_with_filename_and_spectrum(self, ii, dtype,
                                                    input_files, input_spectra):
        fname = input_files[ii]
        spec = input_spectra if issubclass(dtype, Table) else input_spectra[0]
        src = Source(filename=fname, spectra=spec)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, dtype)
        src.shift(0.1, 0.2)

    def test_initialised_with_old_style_arrays(self):
        x, y = [0, 1], [0, -1]
        ref, weight = [0, 0], [1, 10]
        lam = np.linspace(0.5, 2.5, 11) * u.um
        spectra = np.ones(11) * PHOTLAM
        src = Source(x=x, y=y, ref=ref, weight=weight, lam=lam, spectra=spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, Table)
        src.shift(0.1, 0.2)


class TestSourceAddition:
    def test_ref_column_always_references_correct_spectrum(self, table_source,
                                                           image_source):
        image_source.append(table_source)
        comb_refs = image_source.fields[1].field["ref"]
        tbl_refs = table_source.fields[0].field["ref"]
        assert all(tbl_refs.data + 1 == comb_refs.data)
        assert image_source.fields[0].header["SPEC_REF"] == 0
        image_source.shift(0.1, 0.2)

    def test_same_as_above_but_reversed(self, table_source, image_source):
        new_source = table_source + image_source
        comb_refs = new_source.fields[0].field["ref"]
        tbl_refs = table_source.fields[0].field["ref"]
        assert all(tbl_refs.data == comb_refs.data)
        assert new_source.fields[1].header["SPEC_REF"] == 3
        new_source.shift(0.1, 0.2)

    def test_imagehdu_with_empty_spec_ref_is_handled(self, table_source,
                                                     image_source):
        image_source.fields[0].header["SPEC_REF"] = ""
        new_source = table_source + image_source
        assert new_source.fields[1].header["SPEC_REF"] == ""

    def test_fits_image_and_array_image_are_added_correctly(self):
        img_src = so._image_source()
        fits_src = so._fits_image_source()

        img_fits_src = img_src + fits_src
        fits_img_src = fits_src + img_src

        assert len(img_src.fields) == 1
        assert len(fits_src.fields) == 1
        assert len(img_fits_src.fields) == 2
        assert len(fits_img_src.fields) == 2
        assert (fits_img_src.fields[0].data == fits_src.fields[0].data).all()
        assert img_fits_src.fields[0] is not img_src.fields[0]

    @pytest.mark.skip(reason="_meta_dicts was removed, find a better way to perform the same check...")
    def test_meta_data_is_passed_on_when_added(self, table_source, image_source):
        table_source.fields[0].meta["hello"] = "world"
        image_source.fields[0].meta["servus"] = "oida"
        new_source = table_source + image_source

        assert len(new_source.fields) == len(new_source._meta_dicts)
        assert new_source._meta_dicts[0]["hello"] == "world"
        assert new_source._meta_dicts[1]["servus"] == "oida"

    @pytest.mark.skip(reason="_meta_dicts was removed, find a better way to perform the same check...")
    def test_empty_source_is_the_additive_identity(self, image_source):
        new_source_1 = Source() + image_source
        assert len(new_source_1.fields) == len(new_source_1._meta_dicts)
        new_source_2 = image_source + Source()
        assert len(new_source_2.fields) == len(new_source_2._meta_dicts)


@pytest.mark.xfail
class TestSourceShift:
    def test_that_it_does_what_it_should(self):
        assert False


@pytest.mark.xfail
class TestSourceRotate:
    def test_that_it_does_what_it_should(self):
        assert False


class TestSpectraListConverter:
    def test_works_for_arrays(self):
        spec = source_utils.convert_to_list_of_spectra(
                np.array([0, 1, 1, 0]), np.array([1, 2, 3, 4]))
        assert isinstance(spec[0], SourceSpectrum)

    def test_works_for_2d_arrays(self):
        spec = source_utils.convert_to_list_of_spectra(
            np.array([[0, 1, 1, 0], [0, 1, 1, 0]]),
            np.array([1, 2, 3, 4]))
        assert all(isinstance(sp, SourceSpectrum) for sp in spec)

    def test_works_for_multiple_1d_arrays(self):
        spec = source_utils.convert_to_list_of_spectra(
            [np.array([0, 1, 1, 0]), np.array([0, 1, 1, 0])],
            np.array([1, 2, 3, 4]))
        assert all(isinstance(sp, SourceSpectrum) for sp in spec)

    def test_throws_for_array_mismatch(self):
        with pytest.raises(TypeError):
            source_utils.convert_to_list_of_spectra(
                np.array([0, 1, 1, 0]), [1, 2, 3, 4])

    def test_throws_for_multiple_array_mismatch(self):
        with pytest.raises(ValueError):
            source_utils.convert_to_list_of_spectra(
                [np.array([0, 1, 1, 0]), [0, 1, 1, 0]],
                [np.array([1, 2, 3, 4]), [1, 2, 3, 4]])


def test_cube_source_field():
    size = 5
    hdu = fits.ImageHDU(data=np.arange(size**3).reshape(3*(size,)))

    hdu.header["CUNIT1"] = "arcsec"
    hdu.header["CUNIT2"] = "arcsec"
    hdu.header["CUNIT3"] = "um"
    hdu.header["CTYPE3"] = "WAVE"
    hdu.header["CRVAL1"] = 0
    hdu.header["CRVAL2"] = 0
    csf = CubeSourceField(hdu)

    np.testing.assert_equal(csf.waveset.value, np.arange(1, 6))
    csf.shift(2, 3)
    assert csf.header["CRVAL1"] == 2
    assert csf.header["CRVAL2"] == 3

    _, ax = plt.subplots()
    csf.plot(ax, "red")


def test_throws_for_invalid_ref():
    with pytest.raises(KeyError):
        # Minimal Table
        tbl = Table(data=[[0, 1]], names=["ref"])
        TableSourceField(tbl, {0: None})
