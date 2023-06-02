"""
Spectral grating efficiencies
"""
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import units as u

from .effects import Effect
from .ter_curves import TERCurve
from ..utils import find_file

class SpectralEfficiency(Effect):
    """
    Applies the grating efficiency (blaze function) for a SpectralTraceList
    """

    def __init__(self, filename, **kwargs):
        super().__init__(**kwargs)
        self.filename = find_file(filename)
        self.efficiencies = self.get_efficiencies_from_file(self.filename)


    def get_efficiencies_from_file(self, fname):
        """Reads effciencies from file, returns a dictionary"""

        hdul = fits.open(fname)
        efficiencies = {}
        for hdu in hdul[2:]:
            name = hdu.header['EXTNAME']
            lam = hdu.data['wavelength'] * u.um   # check units explicitely
            trans = hdu.data['transmission']
            effic = TERCurve(wavelength=lam, transmission=trans)
            efficiencies[name] = effic

        hdul.close()
        return efficiencies

    def plot(self):
        """Plot the grating efficiencies"""
        for name, effic in self.efficiencies.items():
            wave = effic.throughput.waveset
            plt.plot(wave.to(u.um), effic.throughput(wave), label=name)

        plt.xlabel("Wavelength [um]")
        plt.ylabel("Grating efficiency")
        plt.title(f"Grating efficiencies from {self.filename}")
        plt.legend()
        plt.show()
