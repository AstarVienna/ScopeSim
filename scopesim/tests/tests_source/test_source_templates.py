import pytest
from pytest import approx

from matplotlib import pyplot as plt

from astropy import units as u
from astropy.table import Table

from scopesim import load_example_optical_train
from scopesim.source import source_templates as src_ts
from scopesim.source.source import Source

PLOTS = False


class TestStar:
    @pytest.mark.parametrize("flux", [5, 5*u.ABmag, 36.31*u.Jy])
    def test_accepts_mag_ABmag_jansky(self, flux):
        src = src_ts.star(flux=flux)
        assert src.fields[0]["weight"] == approx(0.01, rel=0.01)
        src.shift(0.1, 0.2)


class TestStarField:
    def test_star_field_return_source_object(self):
        src = src_ts.star_field(100, 15, 25, 60)
        assert isinstance(src, Source)
        src.shift(0.1, 0.2)

    def test_star_field_throws_error_with_no_kwargs(self):
        with pytest.raises(TypeError):
            src_ts.star_field()

    def test_star_fields_data(self):
        src = src_ts.star_field(100, 15, 25, 60)
        assert isinstance(src.fields[0].field, Table)
        assert all(src.fields[0]["weight"] == 10**(-0.4 * src.fields[0]["mag"]))
        src.shift(0.1, 0.2)

    def test_makes_grid_for_jansky_flux(self):
        src = src_ts.star_field(4, 3631*u.Jy, 36.31*u.Jy, 1)
        assert src.fields[0]["weight"][0] == approx(1, rel=0.01)
        assert src.fields[0]["weight"][-1] == approx(0.01, rel=0.01)
        src.shift(0.1, 0.2)


def test_all_zero_spectra_line_up():
    mag = 0
    vega = src_ts.vega_spectrum(mag)
    ab = src_ts.ab_spectrum(mag)
    st = src_ts.st_spectrum(mag)

    wave = 0.55 * u.um
    assert st(wave).value == approx(vega(wave).value, rel=0.03)
    assert ab(wave).value == approx(vega(wave).value, rel=0.03)


class TestUniformIllumination:
    @pytest.mark.usefixtures("protect_currsys")
    def test_makes_source_and_runs_through_basic_instrument(self):
        opt = load_example_optical_train()

        src = src_ts.uniform_illumination(xs=[-50, 50], ys=[-20, 30],
                                          pixel_scale=1, flux=1*u.mag)
        opt.observe(src)
        im = opt.image_planes[0].data

        if PLOTS:
            plt.imshow(im)
            plt.show()

        assert im[512, 512] > 10 * im[0, 0]
        src.shift(0.1, 0.2)

    def test_loads_for_micado_15arcsec_slit(self):
        illum = src_ts.uniform_illumination(xs=[-8, 8], ys=[-0.03, 0.03],
                                            pixel_scale=0.004, flux=1*u.mJy)
        assert isinstance(illum, Source)
        illum.shift(0.1, 0.2)
