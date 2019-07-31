import synphot
from synphot import units, SourceSpectrum, SpectralElement, Observation
from synphot.models import BlackBodyNorm1D, Box1D, Empirical1D
from scopesim.source.source_utils import get_vega_spectrum, rebin_spectra
import numpy as np
#from scopesim import utils
import astropy.units as u
import pytest


vega = get_vega_spectrum()
bb = SourceSpectrum(BlackBodyNorm1D, temperature=5000)

wave = np.arange(500, 10500, 500) * u.Angstrom
flux = np.random.uniform(size=len(wave)) * 1e-17 * units.FLAM
sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, keep_neg=False)

bandpass = SpectralElement(Box1D, amplitude=1, x_0=5500, width=500
                           )
w1 = np.arange(1000, 10000, 0.1)   # linear scale
w2 = np.geomspace(1000, 10000, 1000)  # log scale
w3 = np.arange(30, 100, 0.1)**2   # a quadratic scale
w4 = np.sort(np.concatenate((w2, w3)))  # a widely changing wave scale
w5 = np.linspace(0.1, 1, 1000) * u.um   # testing units

@pytest.mark.parametrize("spectra", [vega, bb, sp])
@pytest.mark.parametrize("new_waveset", [w1, w2, w3, w4, w5])
def test_rebin_spectra(spectra, new_waveset):

    obs1 = Observation(spectra, bandpass)

    new_spec = rebin_spectra(spectra, new_waveset)
    obs2 = Observation(new_spec, bandpass)

    counts1 = obs1.countrate(area=1).value
    counts2 = obs2.countrate(area=1).value

    assert np.isclose(counts1, counts2, rtol=1e-3)
