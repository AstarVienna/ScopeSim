
import synphot
#from synphot import units, SourceSpectrum, SpectralElement, Observation
#from synphot.models import BlackBodyNorm1D, Box1D
from scopesim.source.source_utils import photons_in_range, new_photons_in_range, get_vega_spectrum
import numpy as np
#from scopesim import utils
#import astropy.units as u
import pytest


# Data for the test
spectrum = get_vega_spectrum()
wmin = np.linspace(0.4, 2.2, 11)
wmax = wmin + np.random.rand(len(wmin))*wmin/2
values = [(float(lmin), float(lmax)) for lmin, lmax in zip(wmin,wmax)]

@pytest.mark.parametrize("area", np.logspace(0, 3, 4))
@pytest.mark.parametrize("wmin,wmax", values)
def test_photons_in_range(area, wmin, wmax):
    """
    Test old and new photon in rage on area and wmin, wmax


    :param area:
    :return:
    """

    counts1 = photons_in_range([spectrum], wave_min=wmin, wave_max=wmax, area=area)
    counts2 = new_photons_in_range(spectrum, wave_min=wmin, wave_max=wmax,area=area)

    print("\n")
    print("-"*50)
    print(" wmin:", wmin,
          " wmax:", wmax,
          " area:", area,
          " counts1:", counts1,
          " counts2:", counts2)
    print("-"*50)

    assert np.isclose(counts1.value, counts2.value, rtol=1e-3)



