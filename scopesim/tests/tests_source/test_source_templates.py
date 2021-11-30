import pytest
from pytest import approx
from scopesim.source import source_templates
from _pytest.python_api import approx
from astropy import units as u
from astropy.table import Table

from scopesim.source import source_templates as src_ts
from scopesim.source.source import Source


class TestStar:
    @pytest.mark.parametrize("flux", [5, 5*u.ABmag, 36.31*u.Jy])
    def test_accepts_mag_ABmag_jansky(self, flux):
        src = src_ts.star(flux=flux)
        assert src.fields[0]["weight"] == approx(0.01, rel=0.01)


class TestStarField:
    def test_star_field_return_source_object(self):
        src = src_ts.star_field(100, 15, 25, 60)
        assert isinstance(src, Source)

    def test_star_field_throws_error_with_no_kwargs(self):
        with pytest.raises(TypeError):
            src_ts.star_field()

    def test_star_fields_data(self):
        src = src_ts.star_field(100, 15, 25, 60)
        assert isinstance(src.fields[0], Table)
        assert all(src.fields[0]["weight"] == 10**(-0.4 * src.fields[0]["mag"]))

    def test_makes_grid_for_jansky_flux(self):
        src = src_ts.star_field(4, 3631*u.Jy, 36.31*u.Jy, 1)
        assert src.fields[0]["weight"][0] == approx(1, rel=0.01)
        assert src.fields[0]["weight"][-1] == approx(0.01, rel=0.01)


def test_all_zero_spectra_line_up():
    mag = 0
    vega = source_templates.vega_spectrum(mag)
    ab = source_templates.ab_spectrum(mag)
    st = source_templates.st_spectrum(mag)

    wave = 0.55 * u.um
    assert st(wave).value == approx(vega(wave).value, rel=0.03)
    assert ab(wave).value == approx(vega(wave).value, rel=0.03)
