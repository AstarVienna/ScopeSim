import pytest
from pytest import approx
import numpy as np
from astropy import units as u
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.optics.fov import FieldOfView
from scopesim.optics.fov_utils import get_cube_waveset


PLOTS = False


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


class TestInit:
    def test_initialises_with_nothing_raise_error(self):
        with pytest.raises(TypeError):
            FieldOfView()

    def test_throws_error_if_no_wcs_in_header(self):
        with pytest.raises(ValueError):
            FieldOfView(fits.Header(), (1, 2) * u.um, area=1*u.m**2)

    def test_initialises_with_header_and_waverange(self):
        hdr = _fov_190_210_um().header
        the_fov = FieldOfView(hdr, (1, 2)*u.um, area=1*u.m**2)
        assert isinstance(the_fov, FieldOfView)


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

        assert fov.fields[0].data.shape == (51, 26)
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

        assert fov.fields[0].shape == (3, 51, 26)

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
        assert fov.fields[1].shape == (51, 26)
        assert len(fov.fields[2]) == 2

        assert len(fov.spectra) == 3
        assert fov.fields[1].header["SPEC_REF"] == 0
        for spec in fov.spectra.values():
            assert spec.waveset[0].value == approx(1.97e4)
            assert spec.waveset[-1].value == approx(2.02e4)     # Angstrom

    # Below are tests from original FieldOfView object

    def test_ignores_fields_outside_fov_boundary(self):
        src = so._combined_source(dx=[200, 200, 200])
        src.fields[0]["x"] += 200

        fov = _fov_197_202_um()
        fov.extract_from(src)

        assert len(fov.fields) == 0

    def test_all_spectra_are_referenced_correctly(self):
        src = so._image_source() + so._cube_source() + so._table_source()
        fov = _fov_190_210_um()
        fov.extract_from(src)
        # check the same spectrum object is referenced by both lists
        assert fov.fields[0].header["SPEC_REF"] == \
               src.fields[0].header["SPEC_REF"]
        assert np.all([fov.fields[2][i]["ref"] == \
                       src.fields[2][i]["ref"] for i in range(4)])

        def test_contains_all_fields_inside_fov(self, basic_fov_header,
                                                cube_source,
                                                image_source, table_source):
            src = image_source + cube_source + table_source
            the_fov = FieldOfView(basic_fov_header, (1, 2) * u.um,
                                  area=1 * u.m ** 2)
            the_fov.extract_from(src)
            assert len(the_fov.fields) == 3
            assert isinstance(the_fov.fields[0], fits.ImageHDU)
            assert isinstance(the_fov.fields[1], fits.ImageHDU)
            assert the_fov.fields[1].header["NAXIS"] == 3
            assert isinstance(the_fov.fields[2], Table)


class TestMakeCube:
    def test_makes_cube_from_table(self):
        src_table = so._table_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_table)

        cube = fov.make_cube_hdu()

        in_sum = 0
        waveset = fov.spectra[0].waveset
        for x, y, ref, weight in src_table.fields[0]:
            flux = src_table.spectra[ref](waveset).to(u.ph/u.s/u.m**2/u.um).value
            in_sum += np.sum(flux) * weight

        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum, rel=0.01)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_imagehdu(self):
        src_image = so._image_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_image)

        cube = fov.make_cube_hdu()

        waveset = np.linspace(1.9, 2.1, np.shape(cube)[0]) * u.um
        spec = fov.spectra[0](waveset).to(u.ph/u.s/u.m**2/u.um).value
        in_sum = np.sum(src_image.fields[0].data) * np.sum(spec)
        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum, rel=0.01)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_other_cube_imagehdu(self):
        import scopesim as sim
        sim.rc.__currsys__["!SIM.spectral.spectral_bin_width"] = 0.01
        src_cube = so._cube_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202_um()
        fov.extract_from(src_cube)

        cube = fov.make_cube_hdu()

        # layer 74 to 77 are extracted by FOV
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

        cube = fov.make_cube_hdu()

        # layer 74 to 77 are extracted by FOV
        bin_widths = np.array([0.01, 0.02, 0.01])[:, None, None] * 1e4      # um -> AA
        in_sum = 2 * np.sum(src_cube.fields[0].data[74:77, :, :] * bin_widths)

        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum, rel=1e-3)

        if PLOTS:
            plt.imshow(cube.data[0, :, :], origin="lower")
            plt.show()

    def test_makes_cube_from_all_types_of_source_object(self):
        src_all = so._table_source() + \
                  so._image_source(dx=-4, dy=-4) + \
                  so._cube_source(weight=1e-8, dx=4)

        fov = _fov_190_210_um()
        fov.extract_from(src_all)

        cube = fov.make_cube_hdu()

        # sum up the expected flux in the output cube
        # bin_width * half width edge bin * PHOTLAM -> SI
        waveset = fov.spectra[0].waveset

        table_sum = 0
        for x, y, ref, weight in src_all.fields[0]:
            flux = src_all.spectra[ref](waveset).value
            table_sum += np.sum(flux) * weight

        ref = src_all.fields[1].header["SPEC_REF"]
        spec = fov.spectra[ref](waveset).value
        image_sum = np.sum(src_all.fields[1].data) * np.sum(spec)

        cube_sum = np.sum(src_all.fields[2].data[70:81, :, :]) * 0.02 * 0.95

        in_sum = table_sum + image_sum + cube_sum
        out_sum = np.sum(cube.data)

        assert out_sum == approx(in_sum, rel=0.02)

        if PLOTS:
            im = cube.data[0, :, :]
            plt.imshow(im, origin="lower", norm=LogNorm(), vmin=1e-8)
            plt.show()

    @pytest.mark.parametrize("src", [so._table_source(),
                                     so._image_source(),
                                     so._cube_source()])
    def test_cube_has_full_wcs(self, src):
        fov = _fov_190_210_um()
        fov.extract_from(src)

        cube = fov.make_cube_hdu()

        assert "CDELT3" in cube.header
        assert "CRVAL3" in cube.header
        assert "CRPIX3" in cube.header
        assert "CUNIT3" in cube.header
        assert "CTYPE3" in cube.header


class TestMakeImage:
    def test_makes_image_from_table(self):
        src_table = so._table_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_table)

        img = fov.make_image_hdu()

        in_sum = 0
        waveset = fov.spectra[0].waveset
        for x, y, ref, weight in src_table.fields[0]:
            flux = src_table.spectra[ref](waveset).to(u.ph/u.s/u.m**2/u.um)
            flux *= 1 * u.m**2 * 0.02 * u.um * 0.9      # 0.9 is to catch the half bins at either end
            in_sum += np.sum(flux).value * weight
        out_sum = np.sum(img.data)

        if PLOTS:
            plt.imshow(img.data, origin="lower")
            plt.show()

        assert out_sum == approx(in_sum, rel=0.02)

    def test_makes_image_from_image(self):
        src_image = so._image_source()  # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_image)

        img = fov.make_image_hdu()

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
        src_cube = so._cube_source()  # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202_um()
        fov.extract_from(src_cube)

        image = fov.make_image_hdu()

        # layer 74 to 77 are extracted by FOV
        in_sum = np.sum(src_cube.fields[0].data[74:77, :, :])
        out_sum = np.sum(image.data)

        if PLOTS:
            plt.imshow(image.data, origin="lower")
            plt.show()

        assert out_sum == approx(in_sum)

    def test_makes_image_from_all_types_of_source_object(self):
        src_all = so._table_source() + \
                  so._image_source(dx=-4, dy=-4) + \
                  so._cube_source(weight=1e-8, dx=4)

        fov = _fov_190_210_um()
        fov.extract_from(src_all)

        image = fov.make_image_hdu()

        # sum up the expected flux in the output cube
        waveset = fov.spectra[0].waveset
        scale_factor = 0.02 * 0.91 * 1e8  # bin_width * half width edge bin * PHOTLAM -> SI

        table_sum = 0
        for x, y, ref, weight in src_all.fields[0]:
            flux = src_all.spectra[ref](waveset).value
            table_sum += np.sum(flux) * weight * scale_factor

        ref = src_all.fields[1].header["SPEC_REF"]
        spec = fov.spectra[ref](waveset).value
        image_sum = np.sum(src_all.fields[1].data) * np.sum(spec) * scale_factor

        cube_sum = np.sum(src_all.fields[2].data[70:81, :, :]) * 0.02 * 0.95

        in_sum = table_sum + image_sum + cube_sum
        out_sum = np.sum(image.data)

        assert out_sum == approx(in_sum, rel=0.01)

        if PLOTS:
            im = image.data
            plt.imshow(im, origin="lower", norm=LogNorm(), vmin=1e-8)
            plt.show()

class TestMakeSpectrum:
    def test_make_spectrum_from_table(self):
        src_table = so._table_source()            # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_table)

        spec = fov.make_spectrum()

        in_sum = np.sum([n * spec(fov.waveset).value
                        for n, spec in zip([3, 1, 1], src_table.spectra)])      # sum of weights [3,1,1]
        out_sum = np.sum(spec(fov.waveset).value)

        assert in_sum == approx(out_sum)

    def test_make_spectrum_from_image(self):
        src_image = so._image_source()  # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_190_210_um()
        fov.extract_from(src_image)

        spec = fov.make_spectrum()

        in_sum = np.sum(src_image.fields[0].data) * \
                 np.sum(src_image.spectra[0](fov.waveset).value)
        out_sum = np.sum(spec(fov.waveset).value)

        assert in_sum == approx(out_sum)

    def test_make_spectrum_from_cube(self):
        src_cube = so._cube_source()  # 10x10" @ 0.2"/pix, [0.5, 2.5]m @ 0.02µm
        fov = _fov_197_202_um()
        fov.extract_from(src_cube)

        spec = fov.make_spectrum()

        in_sum = np.sum(src_cube.fields[0].data[74:77, :, :]) * 1e-8
        out_sum = np.sum(spec(fov.waveset).value)

        assert in_sum == approx(out_sum)

    def test_makes_spectrum_from_all_types_of_source_object(self):
        src_table = so._table_source()
        src_image = so._image_source(dx=-4, dy=-4)
        src_cube = so._cube_source(weight=1e-8, dx=4)
        src_all = src_table + src_image + src_cube

        fov = _fov_190_210_um()
        fov.extract_from(src_all)

        spec = fov.make_spectrum()

        table_sum = np.sum([n * spec(fov.waveset).value
                            for n, spec in zip([3, 1, 1], src_table.spectra)])  # sum of weights [3,1,1]
        image_sum = np.sum(src_image.fields[0].data) * \
                    np.sum(src_image.spectra[0](fov.waveset).value)
        cube_sum = np.sum(src_cube.fields[0].data[70:81, :, :]) * 1e-8

        in_sum = table_sum + image_sum + cube_sum
        out_sum = np.sum(spec(fov.waveset).value)

        assert in_sum == approx(out_sum)

        if PLOTS:
            waves = fov.waveset
            plt.plot(waves, spec(waves))
            for spectrum in src_all.spectra:
                plt.plot(waves, spectrum(waves))
            plt.show()

class TestMakeSpectrumImageCubeAllPlayNicely:
    def test_make_cube_and_make_spectrum_return_the_same_fluxes(self):
        src_all = so._table_source() + \
                  so._image_source(dx=-4, dy=-4) + \
                  so._cube_source(weight=1e-8, dx=4)

        fov = _fov_190_210_um()
        fov.extract_from(src_all)

        # if photlam, units of ph / s / cm2 / AA, else units of ph / s / voxel
        cube = fov.make_cube_hdu()
        cube_waves = get_cube_waveset(cube.header)
        cube_spectrum = cube.data.sum(axis=2).sum(axis=1)

        # always units of ph / s / cm-2 / AA-1
        waves = fov.waveset
        spectrum = fov.make_spectrum()(waves).value

        bin_width = 0.02        # um
        photlam_to_si = 1e8     # cm-2 AA-1 --> m-2 um-1
        edge_compensation = 0.91

        spectrum *= bin_width * photlam_to_si * edge_compensation

        if PLOTS:
            plt.plot(waves, spectrum, "k")
            plt.plot(cube_waves, cube_spectrum, "r")
            plt.show()

        assert np.sum(cube.data) == approx(np.sum(spectrum), rel=0.001)

    def test_make_cube_and_make_image_return_the_same_fluxes(self):
        src_all = so._cube_source(weight=1e-8, dx=4)
                  # so._image_source(dx=-4, dy=-4) #
                  # so._table_source()# + \

        fov = _fov_190_210_um()
        fov.extract_from(src_all)

        # if photlam, units of ph / s / cm2 / AA, else units of ph / s / voxel
        cube_sum = np.sum(fov.make_cube_hdu().data)
        # if photlam, units of ph / s / cm2 / AA, else units of ph / s / pixel
        image_sum = np.sum(fov.make_image_hdu().data)

        assert cube_sum == approx(image_sum, rel=0.05)
