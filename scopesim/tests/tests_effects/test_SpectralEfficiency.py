"""Tests for class SpectralEfficiency"""
import os
import pytest

from astropy.io import fits

from scopesim import rc
from scopesim.effects import SpectralEfficiency, TERCurve
from scopesim.utils import find_file

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
if FILES_PATH not in rc.__search_path__:
    rc.__search_path__ += [FILES_PATH]

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

@pytest.fixture(name="speceff", scope="class")
def fixture_speceff():
    """Instantiate SpectralEfficiency object"""
    return SpectralEfficiency(filename="TER_grating.fits")

class TestSpectralEfficiency:
    def test_initialises_from_file(self, speceff):
        assert isinstance(speceff, SpectralEfficiency)

    def test_initialises_from_hdulist(self):
        fitsfile = find_file("TER_grating.fits")
        hdul = fits.open(fitsfile)
        speceff = SpectralEfficiency(hdulist=hdul)
        assert isinstance(speceff, SpectralEfficiency)

    def test_has_efficiencies(self, speceff):
        efficiencies = speceff.efficiencies
        assert all(isinstance(effic, TERCurve)
                   for _, effic in efficiencies.items())
