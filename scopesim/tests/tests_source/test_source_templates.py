import pytest
import scopesim.source.spectrum_templates
from _pytest.python_api import approx
from astropy import units as u
from astropy.table import Table

from scopesim.source import source_templates as src_ts
from scopesim.source.source import Source


class TestStarField:
    def test_star_field_return_source_object(self):
        src = src_ts.star_field(100, 15, 25, 60)
        assert isinstance(src, Source)

    def test_star_field_throws_error_with_no_kwargs(self):
        with pytest.raises(TypeError):
            src_ts.star_field()

    def test_star_fiels_data(self):
        src = src_ts.star_field(100, 15, 25, 60)
        assert isinstance(src.fields[0], Table)
        assert all(src.fields[0]["weight"] == 10**(-0.4 * src.fields[0]["mag"]))


def test_all_zero_spectra_line_up():
    mag = 0
    vega = scopesim.source.spectrum_templates.vega_spectrum(mag)
    ab = scopesim.source.spectrum_templates.ab_spectrum(mag)
    st = scopesim.source.spectrum_templates.st_spectrum(mag)

    wave = 0.55 * u.um
    assert st(wave).value == approx(vega(wave).value, rel=0.03)
    assert ab(wave).value == approx(vega(wave).value, rel=0.03)