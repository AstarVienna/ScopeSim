"""
Spectral grating efficiencies
"""
import logging
import numpy as np
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

from .effects import Effect
from .ter_curves import TERCurve
from .ter_curves_utils import apply_throughput_to_cube
from ..utils import find_file
from ..base_classes import FieldOfViewBase, FOVSetupBase

class SpectralEfficiency(Effect):
    """
    Applies the grating efficiency (blaze function) for a SpectralTraceList
    """

    def __init__(self, filename, **kwargs):
        super().__init__(**kwargs)

        params = {"z_order": [630]}
        self.meta.update(params)

        self.filename = find_file(filename)
        self.efficiencies = self.get_efficiencies_from_file(self.filename)
        print("Hello, this is SpectralEfficiency init")

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


    def apply_to(self, obj, **kwargs):
        """
        Interface between FieldOfView and SpectralEfficiency

        """
        print("Hello, this is SpectralEfficiency.apply_to")
        print(obj.meta['trace_id'])

        if isinstance(obj, FOVSetupBase):
            # I don't think this is needed for the Efficiency - we should get a fully formed FOV
            print("Got FOVSetupBase")
        if isinstance(obj, FieldOfViewBase):
            # Application to field of view
            if obj.cube is None:
                print("Efficiency: no cube")
            if obj.hdu is None:
                print("Efficiency: no hdu")
            else:
                print("Efficiency: hdu", obj.hdu.data.shape)
            trace_id = obj.meta['trace_id']
            try:
                effic = self.efficiencies[trace_id]
            except KeyError:
                logging.warning("No grating efficiency for trace %s" % trace_id)
                return obj
            wcs = WCS(obj.hdu.header).spectral
            wave_cube = wcs.all_pix2world(np.arange(obj.hdu.data.shape[0]), 0)[0]
            wave_cube = (wave_cube * u.Unit(wcs.wcs.cunit[0])).to(u.AA)
            print(wave_cube)
            print(effic.throughput(wave_cube))
            np.savetxt(f"efficcurve_{trace_id}.txt", (wave_cube, effic.throughput(wave_cube)))
            obj.hdu.writeto(f"before_{trace_id}.fits")
            obj.hdu = apply_throughput_to_cube(obj.hdu, effic.throughput)
            obj.hdu.writeto(f"after_{trace_id}.fits")
        return obj

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
