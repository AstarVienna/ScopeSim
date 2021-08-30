import pytest
from pytest import approx
import numpy as np
from astropy import units as u

from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.optics.fov2 import FieldOfView


def _fov_190_210_um():
    hdr = ho._fov_header()  # 20x20" @ 0.2" --> [-10, 10]"
    wav = [1.9, 2.1] * u.um
    fov = FieldOfView(hdr, wav)
    return fov


def _fov_197_202():
    hdr = ho._fov_header()  # 20x20" @ 0.2" --> [-10, 10]"
    wav = [1.97000000001, 2.02] * u.um  # Needs [1.98, 2.00, 2.02] µm --> 3 slices
    fov = FieldOfView(hdr, wav)
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
        src = so._image_source(dx=10)       # 10x10" @ 0.2"/pix
        fov = _fov_190_210_um()
        fov.extract_from(src)

        assert fov.fields[0].data.shape == (51, 25)
        assert len(fov.spectra[0].waveset) == 11
        assert fov.spectra[0].waveset[0].value == approx(19000)

    def test_extract_3d_cube_from_hduimage(self):
        src = so._cube_source()             # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202()
        fov.extract_from(src)

        s198, s200, s202 = fov.fields[0].data.sum(axis=(2,1))
        assert s198 == approx(s200, rel=0.02)
        assert s202 == approx(s200 * 0.5, rel=0.02)
        assert fov.fields[0].data.shape == (3, 51, 51)
        assert len(fov.spectra) == 0

    def test_extract_3d_cube_that_is_offset_relative_to_fov(self):
        src = so._cube_source(dx=10)        # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm, centre offset to (10, 0)"
        fov = _fov_197_202()
        fov.extract_from(src)

        assert fov.fields[0].shape == (3, 51, 25)


    def test_extract_one_of_each_type_from_source_object(self):
        src_table = so._table_source()              # 4 sources, put two outside of FOV
        src_table.fields[0]["x"] = [-15,-5,0,0] * u.arcsec
        src_table.fields[0]["y"] = [0,0,5,15] * u.arcsec
        src_image = so._image_source(dx=10)         # 10x10" @ 0.2"/pix
        src_cube = so._cube_source()                # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        src = src_cube + src_image + src_table

        fov = _fov_197_202()
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
        pass

    def test_makes_cube_from_imagehdu(self):
        pass

    def test_makes_cube_from_other_cube_imagehdu(self):
        pass
