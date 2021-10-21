import pytest
from pytest import approx
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.optics.fov2 import FieldOfView


PLOTS = True


def _fov_190_210_um():
    """ A FOV compatible with 11 slices of so._cube_source()"""
    hdr = ho._fov_header()  # 20x20" @ 0.2" --> [-10, 10]"
    wav = [1.9, 2.1] * u.um
    fov = FieldOfView(hdr, wav, area=1 * u.m ** 2)
    return fov


def _fov_197_202_um():
    """ A FOV compatible with 3 slices of so._cube_source()"""
    hdr = ho._fov_header()  # 20x20" @ 0.2" --> [-10, 10]"
    wav = [1.97000000001, 2.02] * u.um  # Needs [1.98, 2.00, 2.02] µm --> 3 slices
    fov = FieldOfView(hdr, wav, area=1*u.m**2)
    return fov


class TestExtractFrom:
    def test_extract_point_sources_from_table(self):
        src = so._table_source()
        src.fields[0]["x"] = [-15,-5,0,0] * u.arcsec
        src.fields[0]["y"] = [0,0,5,15] * u.arcsec
        fov = _fov_190_210_um()
        fov.extract_from(src)

        assert len(fov.fields[0]) == 2
        assert len(fov.spectra[0].waveset) == 11
        assert fov.spectra[0].waveset[0].value == approx(19000)

    def test_extract_2d_image_from_hduimage(self):
        src = so._image_source(dx=10)       # 10x10" @ 0.2"/pix, offset by 10"
        fov = _fov_190_210_um()
        fov.extract_from(src)

        assert fov.fields[0].data.shape == (51, 25)
        assert len(fov.spectra[0].waveset) == 11
        assert fov.spectra[0].waveset[0].value == approx(19000)

    def test_extract_3d_cube_from_hduimage(self):
        src = so._cube_source()             # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202_um()
        fov.extract_from(src)

        s198, s200, s202 = fov.fields[0].data.sum(axis=(2,1))
        assert s198 == approx(s200, rel=0.02)
        assert s202 == approx(s200 * 0.5, rel=0.02)
        assert fov.fields[0].data.shape == (3, 51, 51)
        assert len(fov.spectra) == 0

    def test_extract_3d_cube_that_is_offset_relative_to_fov(self):
        src = so._cube_source(dx=10)        # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm, centre offset to (10, 0)"
        fov = _fov_197_202_um()
        fov.extract_from(src)

        assert fov.fields[0].shape == (3, 51, 25)

    def test_extract_one_of_each_type_from_source_object(self):
        src_table = so._table_source()              # 4 sources, put two outside of FOV
        src_table.fields[0]["x"] = [-15,-5,0,0] * u.arcsec
        src_table.fields[0]["y"] = [0,0,5,15] * u.arcsec
        src_image = so._image_source(dx=10)         # 10x10" @ 0.2"/pix
        src_cube = so._cube_source()                # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        src = src_cube + src_image + src_table

        fov = _fov_197_202_um()
        fov.extract_from(src)

        assert fov.fields[0].shape == (3, 51, 51)
        assert fov.fields[1].shape == (51, 25)
        assert len(fov.fields[2]) == 2

        assert len(fov.spectra) == 3
        assert fov.fields[1].header["SPEC_REF"] == 0
        for spec in fov.spectra.values():
            assert spec.waveset[0].value == approx(1.97e4)
            assert spec.waveset[-1].value == approx(2.02e4)     # Angstrom


class TestMakeCube:
    def test_makes_cube_from_table(self):
        src_table = so._table_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_table)

        cube = fov.make_cube()

        in_sum = 0
        waveset = fov.spectra[0].waveset
        for x, y, ref, weight in src_table.fields[0]:
            flux = src_table.spectra[ref](waveset).value
            in_sum += np.sum(flux) * weight
        out_sum = np.sum(cube.data)
        assert out_sum == approx(in_sum)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_imagehdu(self):
        src_image = so._image_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_image)

        cube = fov.make_cube()

        waveset = np.linspace(1.9, 2.1, np.shape(cube)[0]) * u.um
        spec = fov.spectra[0](waveset).value
        in_sum = np.sum(src_image.fields[0].data) * np.sum(spec)
        out_sum = np.sum(cube.data)
        assert out_sum == approx(in_sum)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_other_cube_imagehdu(self):
        src_cube = so._cube_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202_um()
        fov.extract_from(src_cube)

        cube = fov.make_cube()

        n = 74      # layer 74 to 77 are extracted by FOV
        in_sum = np.sum(src_cube.fields[0].data[74:77, :, :])
        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_two_similar_cube_imagehdus(self):
        src_cube = so._cube_source() + so._cube_source(dx=1)            # 2 cubes 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202_um()
        fov.extract_from(src_cube)

        cube = fov.make_cube()

        n = 74      # layer 74 to 77 are extracted by FOV
        in_sum = 2 * np.sum(src_cube.fields[0].data[74:77, :, :])
        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_all_types_of_source_object(self):
        src_all = so._cube_source(weight=1e-8, dx=4) + \
                  so._image_source(dx=-4, dy=-4) + \
                  so._table_source()
        fov = _fov_190_210_um()
        fov.extract_from(src_all)

        cube = fov.make_cube()

        # sum up the expected flux in the output cube
        table_sum = 0
        waveset = fov.spectra[0].waveset
        for x, y, ref, weight in src_all.fields[2]:
            flux = src_all.spectra[ref](waveset).value
            table_sum += np.sum(flux) * weight

        ref = src_all.fields[1].header["SPEC_REF"]
        spec = fov.spectra[ref](waveset).value
        image_sum = np.sum(src_all.fields[1].data) * np.sum(spec)

        cube_sum = np.sum(src_all.fields[0].data[70:81, :, :])

        in_sum = table_sum + image_sum + cube_sum
        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower", norm=LogNorm())
            plt.show()


class TestMakeImage:
    def test_makes_image_from_table(self):
        src_table = so._table_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_table)

        img = fov.make_image()

        in_sum = 0
        waveset = fov.spectra[0].waveset
        for x, y, ref, weight in src_table.fields[0]:
            flux = src_table.spectra[ref](waveset).value
            in_sum += np.sum(flux) * weight
        out_sum = np.sum(img.data)
        assert out_sum == approx(in_sum)

        if PLOTS:
            plt.imshow(img.data, origin="lower")
            plt.show()

    def test_makes_image_from_image(self):
        src_image = so._image_source()  # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_image)

        img = fov.make_image()

        src_im_sum = np.sum(src_image.fields[0].data)
        src_spec = src_image.spectra[0](fov.waveset).to(u.ph/u.s/u.m**2/u.um)
        src_flux = 0.91 * np.sum(src_spec * 1 * u.m**2 * 0.02 * u.um).value
        in_sum = src_im_sum * src_flux
        out_sum = np.sum(img.data)

        # 2% discrepancy is because the edge bins in FOV (not of equal value)
        # are multiplied by 0.5 when summed. This test simply multiplies by 0.91
        # to account for the last (11th) bin edge in the range 1.9 : 2.1 : 0.01

        assert out_sum == approx(in_sum, rel=0.01)

        if PLOTS:
            plt.imshow(img.data, origin="lower")
            plt.show()

    def test_makes_image_from_cube(self):
