"""Classes for the METIS Warm Calibration Unit"""

import numpy as np
from astropy.table import Table
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii as ioascii
from scipy.interpolate import interp1d
import yaml

from ..ter_curves import TERCurve
from ...utils import get_logger, seq, find_file,\
    convert_table_comments_to_dict
from ...source.source import Source

logger = get_logger(__name__)

class BlackBodySource(TERCurve):
    """Black Body Source

    This class returns a TERCurve that describes the output emission
    from the integrating sphere of the METIS WCU. The calculations include
    - black-body emission from the black-body source
    - coupling into the integrating sphere
    - integrating sphere magnification factor
    - thermal emission from the integrating sphere
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "z_order": [113, 513],
            "action": "emissivity",
            "position": 0,  # position in surface table
        }
        self.meta.update(params)
        self.meta.update(kwargs)
        if 'config_file' in self.meta:
            with open(find_file(self.meta['config_file'])) as fd:
                config = yaml.safe_load(fd)
                self.meta.update(config)

        # Check on the presence of one vital parameter
        if "rho_tube" not in self.meta:
            raise ValueError("Parameters not present: please provide config file or parameter values"
                             "for BlackBodySource")

        self._background_source = None

        # Load components for the source
        self.bb_to_is = self.bb_to_is_throughput()
        self.rho_tube = get_reflectivity(self.meta['rho_tube'])
        self.rho_is = get_reflectivity(self.meta['rho_is'])
        self.compute_emission()

    @property
    def emission(self):
        return self.surface.emission

    @property
    def background_source(self):
        """Define a source field for the FP mask"""
        if self._background_source is None:
            flux = self.emission
            bg_hdu = fits.ImageHDU()

            hdr = {"BG_SRC": True,
                   "BG_SURF": self.display_name,
                   "CTYPE1": "LINEAR",
                   "CTYPE2": "LINEAR",
                   "CRPIX1": 1024.5,
                   "CRPIX2": 1024.5,
                   "CRVAL1": 0.,
                   "CRVAL2": 0.,
                   "CUNIT1": "ARCSEC",
                   "CUNIT2": "ARCSEC",
                   "CDELT1": 0.00547,
                   "CDELT2": 0.00547,
                   #"BUNIT": "PHOTLAM arcsec-2",
                   "SOLIDANG": "arcsec-2"}
            bg_hdu.header.update(hdr)
            bg_hdu.data = np.zeros((2048, 2048))
            bg_hdu.data[1024, 1024] = 1
            bg_hdu.data[512, 512] = 1
            self._background_source = [Source(image_hdu=bg_hdu, spectra=flux)]

            bg2_hdu = fits.ImageHDU()
            bg2_hdu.header.update(hdr)
            bg2_hdu.data = np.ones((2048, 2048))
            bg2_hdu.data[1024, 1024] = 0
            bg2_hdu.data[512, 512] = 0
            self._background_source.append(Source(image_hdu=bg2_hdu, spectra=0.02 * flux))

        return self._background_source


    def set_temperature(self, bb_temp: [float | u.Quantity]=None,
                        wcu_temp: [float | u.Quantity]=None,
                        is_temp: [float | u.Quantity]=None):
        """Change the black-body temperature

        Parameters
        ----------
        bb_temp, wcu_temp : float, Quantity
            new temperatures for the BB source and the ambient WCU, respectively.
            If float, the unit is assumed to be Kelvin.
        """
        if bb_temp is not None:
            if isinstance(bb_temp, (int, float)):
                bb_temp = bb_temp << u.K

            bb_temp = bb_temp.to(u.K, equivalencies=u.temperature())
            if bb_temp >= 0:
                self.meta["bb_temp"] = bb_temp
            else:
                raise ValueError("bb_temp below absolute zero, not changed")

        if wcu_temp is not None:
            if isinstance(wcu_temp, (int, float)):
                wcu_temp = wcu_temp << u.K

            wcu_temp = wcu_temp.to(u.K, equivalencies=u.temperature())
            if wcu_temp >= 0:
                self.meta["wcu_temp"] = wcu_temp
            else:
                raise ValueError("wcu_temp below absolute zero, not changed")

        if is_temp is not None:
            if isinstance(is_temp, (int, float)):
                is_temp = is_temp << u.K

            is_temp = is_temp.to(u.K, equivalencies=u.temperature())
            if is_temp >= 0:
                self.meta["is_temp"] = is_temp
            else:
                raise ValueError("is_temp below absolute zero, not changed")

        self.compute_emission()


    def bb_to_is_throughput(self):
        """Load throughput and return interpolation function

        This loads a lookup table for the transmission from the black-body
        source to the entrance port of the integrating sphere. It returns
        an interpolation function that can be evaluated on a reflectivity for the tube.

        The model used is the general model, currently with no gaps.
        """

        path = find_file(self.meta['bb_to_is'])
        if path is None:
            return lambda x: 1
        hdul = fits.open(path)
        rho_tube = hdul[1].data['rho_tube']
        throughput = hdul[1].data['t_gen_no_gap']    # TODO: update
        return interp1d(rho_tube, throughput)


    def compute_emission(self):
        """Compute the emission at the exit of the integrating sphere"""
        self.d_is = self.meta["diam_is"] << u.mm
        self.d_is_in = self.meta["diam_is_in"] << u.mm
        self.d_is_out = self.meta["diam_is_out"] << u.mm
        self.bb_temp = self.meta["bb_temp"] << u.K
        self.is_temp = self.meta["is_temp"] << u.K
        self.wcu_temp = self.meta["wcu_temp"] << u.K
        self.emiss_bb = self.meta["emiss_bb"]

        tbl = Table()
        lam = seq(2.2, 15, 0.01) * u.um

        mult_is = self.is_multiplication(lam)

        # continuum black-body source
        bb_scale = 1 * u.ph / (u.s * u.m**2 * u.sr * u.um)
        self.is_bb = BlackBody(self.bb_temp, scale=bb_scale)
        self.flux_bb = (self.emiss_bb * self.is_bb(lam)
                        * (np.pi * self.d_is_in**2 / 4) * (np.pi * u.sr))
        self.flux_bb *= self.bb_to_is(self.rho_tube(lam))
        self.intens_bb = self.flux_bb / (np.pi * self.d_is**2) * mult_is / (np.pi * u.sr)

        # background emission from integrating sphere
        self.is_bg = BlackBody(self.is_temp, scale=bb_scale)
        #self.flux_bg = ((1 - self.rho_is(lam)) * self.is_bg(lam)
        #                * (np.pi * self.d_is**2) * (np.pi * u.sr))
        #self.intens_bg = self.flux_bg / (np.pi * self.d_is**2) * mult_is / (np.pi * u.sr)
        self.intens_bg  = self.rho_is(lam) * self.is_bg(lam)

        self.intensity = self.intens_bb + self.intens_bg
        self.wavelength = lam
        tbl.add_column(lam, name="wavelength")
        tbl.add_column(np.ones_like(lam).value, name="transmission")
        tbl.add_column(self.intensity, name="emission")
        tbl.meta["wavelength_unit"] = tbl["wavelength"].unit
        tbl.meta["emission_unit"] = tbl["emission"].unit

        self.surface.table = tbl
        self.surface.meta.update(tbl.meta)


    def is_multiplication(self, wavelength, nport=2):
        """Intensity multiplication factor for the integrating sphere

        Parameters
        ----------
        wavelength : nd.array, float
            Wavelengths [mm] on which to compute the factor.
        nport: int
            Number of ports that contribute to "missing" reflective area
            on the inside of the integrating sphere. Typical values for
            the METIS WCU are 2 or 4, depending on whether the extra ports
            are counted. Value 0 can be used for testing purposes.
        """
        rho = self.rho_is(wavelength)
        area_in = (self.d_is_in**2
                   * ( 1 - np.cos(np.arcsin(self.d_is_in / self.d_is))))
        area_out = (self.d_is_out**2
                    * (1 - np.cos(np.arcsin(self.d_is_out / self.d_is))))
        area_is = np.pi * self.d_is**2
        if nport == 0:
            frac = ((area_in + (nport - 1) * area_out) / area_is).decompose().value
        else:
            frac = 0

        return rho / (1 - rho * (1 - frac))

    def __str__(self) -> str:
        return f"""{self.__class__.__name__}: "{self.display_name}"
        BlackBody temperature:   {self.meta['bb_temp']}
        Integrating sphere temp: {self.meta['is_temp']}
        WCU temperature:         {self.meta['wcu_temp']}"""


# TODO: put into metis_wcu_utils.py
def get_reflectivity(file_or_number):
    """
    Get a reflectivity from either a file or a number

    This returns a function of wavelength.
    Does not work with Quantity, assumes wavelength in um.
    """
    if isinstance(file_or_number, str):
        tbl = ioascii.read(find_file(file_or_number))
        return interp1d(tbl['wavelength'], tbl['reflection'])

    return lambda x: file_or_number * np.ones_like(x.value)
