"""Tests for class SpectralEfficiency"""

import pytest

from astropy.io import fits

from scopesim.effects import SpectralEfficiency, TERCurve


@pytest.fixture(name="speceff", scope="class")
def fixture_speceff(mock_path):
    """Instantiate SpectralEfficiency object"""
    return SpectralEfficiency(filename=str(mock_path / "TER_grating.fits"))


class TestSpectralEfficiency:
    def test_initialises_from_file(self, speceff):
        assert isinstance(speceff, SpectralEfficiency)

    def test_initialises_from_hdulist(self, mock_path):
        # fitsfile = find_file("TER_grating.fits")
        fitsfile = mock_path / "TER_grating.fits"
        with fits.open(fitsfile) as hdul:
            speceff = SpectralEfficiency(hdulist=hdul)
        assert isinstance(speceff, SpectralEfficiency)

    def test_has_efficiencies(self, speceff):
        efficiencies = speceff.efficiencies
        assert all(isinstance(effic, TERCurve)
                   for effic in efficiencies.values())
