"""Classes for the METIS Warm Calibration Unit"""

import numpy as np
from astropy.table import Table
from astropy.modeling.models import BlackBody, Gaussian1D
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.io import ascii as ioascii
from scipy.interpolate import interp1d
import yaml

from ..ter_curves import TERCurve, FilterCurve
from ...utils import get_logger, seq, find_file,\
    from_currsys
from ...source.source import Source
from ...optics.surface import SpectralSurface


logger = get_logger(__name__)

class WCUSource(TERCurve):
    """Warm Calibration Unit Source

    This class returns a TERCurve that describes the output emission
    from the integrating sphere of the METIS WCU. The calculations include
    - black-body emission from the black-body source
    - coupling into the integrating sphere
    - integrating sphere magnification factor
    - thermal emission from the integrating sphere

    Configuration file
    ------------------
    The default configuration file for the METIS WCU is provided in the irdb.
    To modify the configuration, copy the file <irdb>/METIS/metis_wcu_config.yaml
    to the working directory and edit it there. To use the file do e.g.
    >>> cmds = sim.UserCommands(use_instrument="METIS", set_modes=["wcu_img_lm"]
    >>> cmds["!WCU.config_file"] = "my_wcu_config.yaml"
    >>> metis = sim.OpticalTrain(cmds)

    Changing important parameters
    -----------------------------
    The temperatures of the black-body source, the integrating sphere as well as the
    ambient temperature of the WCU (which is the temperature of the focal-plane mask)
    can be set by calling
    >>> metis['bb_source'].set_temperature(bb_temp=1200*u.K, is_temp=320*u.K,
                                           wcu_temp=295*u.K)

    The focal-plane mask can be changed by calling
    >>> metis['bb_source'].set_mask("pinhole")
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
            config_file = from_currsys(self.meta['config_file'], self.cmds)
            with open(find_file(config_file)) as fd:
                config = yaml.safe_load(fd)
                self.meta.update(config)

        # Check on the presence of one vital parameter
        if "rho_tube" not in self.meta:
            raise ValueError("Parameters not present: please provide config file or parameter values"
                             "for WCUSource")

        self._background_source = None

        # Define surfaces (the IS surface is currently taken from parent class)
        # we add another one here (TODO: make this nicer)
        self.mask_surf = SpectralSurface(cmds = self.cmds)
        self.mask_surf.meta.update(self.meta)

        # Load components for the source
        #self.wavelength = seq(2.2, 15, 0.01) * u.um   # TODO cleverer
        self.get_wavelength()
        self.bb_scale = 1 * u.ph / (u.s * u.m**2 * u.arcsec**2 * u.um)
        self.bb_to_is = self.bb_to_is_throughput()
        self.rho_tube = get_reflectivity(self.meta['rho_tube'])
        self.rho_is = get_reflectivity(self.meta['rho_is'])
        self.rho_mask = get_reflectivity(self.meta['rho_mask'])
        self.compute_lamp_emission()
        self.compute_fp_emission()

    @property
    def emission(self):
        return self.surface.emission

    @property
    def mask_emission(self):
        return self.mask_surf.emission

    @property
    def background_source(self):
        """Define a source field for the FP mask"""

        bb_flux = self.emission
        mask_flux = self.mask_emission
        self._background_source = []
        hdr = {"BG_SRC": True,
               "BG_SURF": self.display_name,
               "CTYPE1": "LINEAR",
               "CTYPE2": "LINEAR",
               "CRPIX1": 1024.5,
               "CRPIX2": 1024.5,
               "CRVAL1": 0.,
               "CRVAL2": 0.,
               "CUNIT1": "arcsec",
               "CUNIT2": "arcsec",
               "CDELT1": 0.00547,
               "CDELT2": 0.00547,
               "BUNIT": "PHOTLAM arcsec-2",
               "SOLIDANG": "arcsec-2"}
        if self.meta['current_mask'] == 'open':
            bg_hdu = fits.ImageHDU()
            bg_hdu.header.update(hdr)
            self._background_source.append(Source(image_hdu=bg_hdu, spectra=bb_flux))
        else:   # TODO: properly define masks
            holearea = 4.45610478e-05 * u.arcsec**2 # LM:  25 um
            # holearea = 0.00031057 * u.arcsec**2   # N:  66 um

            bg_hdu = fits.ImageHDU()
            bg_hdu.header.update(hdr)
            bg_hdu.data = np.zeros((2048, 2048))
            bg_hdu.data[1024, 1024] = 1
            bg_hdu.data[512, 512] = 1
            # When spectrum is intensity (../arcsec2), the image must contain the pixel
            # area (arcsec2). For a true backgroud field, this is taken care of by
            # fov._calc_area_factor, but this is not applied to image data, so we
            # have to do it here.
            bg_hdu.data = bg_hdu.data * holearea  # Should actually be the size of the hole
            self._background_source.append(Source(image_hdu=bg_hdu, spectra=bb_flux))

            pixarea = (hdr['CDELT1'] * u.Unit(hdr['CUNIT1'])
                       * hdr['CDELT2'] * u.Unit(hdr['CUNIT2']))

            bg2_hdu = fits.ImageHDU()
            bg2_hdu.header.update(hdr)
            bg2_hdu.data = np.ones((2048, 2048))
            bg2_hdu.data[1024, 1024] = 0
            bg2_hdu.data[512, 512] = 0
            bg2_hdu.data = bg2_hdu.data * pixarea
            self._background_source.append(Source(image_hdu=bg2_hdu, spectra=mask_flux))
            #self._background_source = [Source(image_hdu=bg2_hdu, spectra=bb_flux)]  # TEST!!!

        return self._background_source

    def set_lamp(self, lamp='bb'):
        """Change the WCU lamp

        The list of available lamps can be obtained with
        >>> metis['wcu_source'].lamps
        The lamp currently in use is
        >>> metis['wcu_source'].current_lamp
        It is changed with, e.g.,
        >>> metis['wcu_source'].set_lamp('laser')
        """
        if lamp not in self.lamps:
            raise ValueError(f"'lamp' needs to be one of {self.lamps}")

        self.meta['current_lamp'] = lamp
        #self.get_wavelength()   # does not depend on the lamp
        self.compute_lamp_emission()
        self.compute_fp_emission()


    def get_wavelength(self):
        """Try to set the appropriate wavelength vector for the mode and filter"""
        if 'wcu_lms' in self.cmds['!OBS.modes']:     ## Need to provide for wcu_lms_extended
            lamc = self.cmds['!OBS.wavelen']
            dlam = self.cmds['!SIM.spectral.spectral_bin_width']
            lam = seq(lamc - 3000 * dlam, lamc + 3000 * dlam, dlam) * u.um
        else:
            filter_name = self.cmds['!OBS.filter_name']
            filename_format = self.cmds['!INST.filter_file_format']
            tempfilter = FilterCurve(filter_name=filter_name,
                                     filename_format=filename_format)
            lammin, lammax = tempfilter.throughput.waverange << u.um
            dlam = self.cmds['!SIM.spectral.spectral_bin_width']
            lam = seq(lammin.value, lammax.value, dlam) * u.um
        self.wavelength = lam


    @property
    def lamps(self):
        """List of available lamps"""
        return self.meta['lamps']

    @property
    def current_lamp(self):
        """Name of the lamp currently in use"""
        return self.meta['current_lamp']

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

        self.compute_lamp_emission()
        self.compute_fp_emission()


    def set_mask(self, fpmask: str):
        """Change the focal-plane mask"""
        masklist = self.meta['fpmasks']
        if fpmask not in masklist:
            raise ValueError(f"fpmask must be one of {masklist}")
        self.meta["current_mask"] = fpmask

        self.compute_lamp_emission()
        self.compute_fp_emission()

    @property
    def current_mask(self):
        """Name of the mask currently in use"""
        return self.meta['current_mask']


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


    def compute_lamp_emission(self):
        """Compute the emission at the exit of the integrating sphere"""
        self.d_is = self.meta["diam_is"] << u.mm
        self.d_is_in = self.meta["diam_is_in"] << u.mm
        self.d_is_out = self.meta["diam_is_out"] << u.mm
        self.bb_temp = self.meta["bb_temp"] << u.K
        self.is_temp = self.meta["is_temp"] << u.K
        self.wcu_temp = self.meta["wcu_temp"] << u.K
        self.emiss_bb = self.meta["emiss_bb"]      # <<<<<< that needs to be a function

        lam = self.wavelength

        mult_is = self.is_multiplication(lam)

        # continuum black-body source
        if self.current_lamp == "bb":
            self.is_lamp = BlackBody(self.bb_temp, scale=self.bb_scale)
            self.flux_lamp = (self.emiss_bb * self.is_lamp(lam)
                            * (np.pi * self.d_is_in**2 / 4) * (np.pi * u.sr))
            self.flux_lamp *= self.bb_to_is(self.rho_tube(lam))
            self.intens_lamp = self.flux_lamp / (np.pi * self.d_is**2) * mult_is / (np.pi * u.sr)
        elif self.current_lamp == "laser":
            self.intens_lamp = self._laser_intensity()
        elif self.current_lamp == "none":
            self.intens_lamp = np.zeros(len(lam)) * self.bb_scale
        else:
            raise ValueError(f"Unknown lamp: {self.current_lamp}")

        # background emission from integrating sphere
        self.is_bg = BlackBody(self.is_temp, scale=self.bb_scale)
        #self.flux_bg = ((1 - self.rho_is(lam)) * self.is_bg(lam)
        #                * (np.pi * self.d_is**2) * (np.pi * u.sr))
        #self.intens_bg = self.flux_bg / (np.pi * self.d_is**2) * mult_is / (np.pi * u.sr)
        self.intens_bg  = self.rho_is(lam) * self.is_bg(lam)

        self.intensity = self.intens_lamp + self.intens_bg

        tbl = Table()
        tbl.add_column(lam, name="wavelength")
        tbl.add_column(np.ones_like(lam).value, name="transmission")
        tbl.add_column(self.intensity, name="emission")
        tbl.meta["wavelength_unit"] = tbl["wavelength"].unit
        tbl.meta["emission_unit"] = tbl["emission"].unit

        self.surface.table = tbl
        self.surface.meta.update(tbl.meta)

    def _laser_intensity(self):
        lam = self.wavelength
        dlam = lam[1] - lam[0]

        power = 5e-3 * u.W / (c.c * c.h / lam) * u.ph
        lamc = 3.39 * u.um
        sigma = 2 * dlam #* u.um
        amp = 1/(sigma * np.sqrt(2 * np.pi))
        line = Gaussian1D(amplitude=amp, mean=lamc, stddev=sigma)
        flux = power * line(lam) / (np.pi * self.d_is_in**2 / 4)
        intens = flux / (np.pi * u.sr)
        return intens.to(self.bb_scale)

    def compute_fp_emission(self):
        """Compute the emission from the opaque part of the focal-plane mask"""
        self.wcu_temp = self.meta["wcu_temp"] << u.K
        self.emiss_mask = 1 - self.meta["rho_mask"]       # <<<<<< that needs to be a function

        lam = self.wavelength

        # continuum black-body source
        #self.bb_scale = 1 * u.ph / (u.s * u.m**2 * u.arcsec**2 * u.um)
        self.mask_em = BlackBody(self.wcu_temp, scale=self.bb_scale)
        self.intens_fp = self.emiss_mask * self.mask_em(lam)

        tbl = Table()
        tbl.add_column(lam, name="wavelength")
        tbl.add_column(np.ones_like(lam).value, name="transmission")
        tbl.add_column(self.intens_fp, name="emission")
        tbl.meta["wavelength_unit"] = tbl["wavelength"].unit
        tbl.meta["emission_unit"] = tbl["emission"].unit

        self.mask_surf.table = tbl
        self.mask_surf.meta.update(tbl.meta)


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
        Current lamp:            {self.current_lamp}     {self.lamps}
        BlackBody temperature:   {self.meta['bb_temp']}
        Integrating sphere temp: {self.meta['is_temp']}
        WCU temperature:         {self.meta['wcu_temp']}
        Focal-plane mask:        {self.meta['current_mask']}  {self.meta['fpmasks']}"""


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
