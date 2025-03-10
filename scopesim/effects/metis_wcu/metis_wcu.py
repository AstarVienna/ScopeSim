"""Classes for the METIS Warm Calibration Unit"""

from typing import ClassVar

import numpy as np
from astropy.table import Table
from astropy.modeling.models import BlackBody, Gaussian1D
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.io import ascii as ioascii
from scipy.interpolate import interp1d
import yaml

from .fpmask import FPMask
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
    >>> metis['wcu_source'].set_temperature(bb_temp=1200*u.K, is_temp=320*u.K,
                                            wcu_temp=295*u.K)

    The flux-controlling masks are currently implemented as a float (in [0, 1]) that
    gives the fraction of flux transmitted. It can be change with
    >>> metis['wcu_source'].set_bb_aperture(0.8)

    The focal-plane mask and/or its orientation and position can be changed by calling
    >>> metis['wcu_source'].set_fpmask("pinhole", angle=10, shift=(0.1, 0))

    To use the lasers instead of the black-body source, do
    >>> metis['wcu_source'].set_lamp("laser")
    Which laser is seen depends on the wavelength of the observation. Note that the
    tunable laser currently cannot be tuned.
    """

    z_order: ClassVar[tuple[int, ...]] = (113, 513)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "action": "emissivity",
            "position": 0,  # position in surface table
        }
        self.meta.update(params)
        self.meta.update(kwargs)
        if 'config_file' in self.meta:
            config_file = from_currsys(self.meta['config_file'], self.cmds)
            with open(find_file(config_file), encoding="utf-8") as fd:
                config = yaml.safe_load(fd)
                self.meta.update(config)

        # Check on the presence of one vital parameter
        if "rho_tube" not in self.meta:
            raise ValueError(
                "Parameters not present: please provide config file or parameter values"
                "for WCUSource")

        self.set_fpmask(maskname=self.meta['current_fpmask'],
                        angle=self.meta["fpmask_angle"],
                        shift=self.meta["fpmask_shift"])
        self._background_source = None

        self.bb_aperture = self.meta['bb_aperture']

        # Define surfaces (the IS surface is currently taken from parent class)
        # we add another one here (TODO: make this nicer)
        self.mask_surf = SpectralSurface(cmds = self.cmds)
        self.mask_surf.meta.update(self.meta)

        # Load components for the source
        self.get_wavelength()
        self.bb_scale = 1 * u.ph / (u.s * u.m**2 * u.arcsec**2 * u.um)
        self.bb_to_is = self.bb_to_is_throughput()
        self.rho_tube = get_reflectivity(self.meta['rho_tube'])
        self.rho_is = get_reflectivity(self.meta['rho_is'])
        self.emiss_mask = self.meta['emiss_mask']

        # Compute the emission components
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

        self._background_source.append(Source(image_hdu=self.fpmask.holehdu,
                                              spectra=bb_flux))
        if self.fpmask.opaquehdu is not None:   # not the open mask
            self._background_source.append(Source(image_hdu=self.fpmask.opaquehdu,
                                                  spectra=mask_flux))
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
                        is_temp: [float | u.Quantity]=None,
                        wcu_temp: [float | u.Quantity]=None):
        """Change the black-body temperature

        Parameters
        ----------
        bb_temp, wcu_temp, is_temp : float, Quantity
            new temperatures for the BB source, the integrating sphere and the
            ambient WCU, arespectively.
            If float, the unit is assumed to be Kelvin.
        """
        if bb_temp is not None:
            with u.set_enabled_equivalencies(u.temperature()):
                bb_temp <<= u.K
            if bb_temp >= 0:
                self.meta["bb_temp"] = bb_temp
            else:
                raise ValueError("bb_temp below absolute zero, not changed")

        if wcu_temp is not None:
            with u.set_enabled_equivalencies(u.temperature()):
                wcu_temp <<= u.K
            if wcu_temp >= 0:
                self.meta["wcu_temp"] = wcu_temp
            else:
                raise ValueError("wcu_temp below absolute zero, not changed")

        if is_temp is not None:
            with u.set_enabled_equivalencies(u.temperature()):
                is_temp <<= u.K
            if is_temp >= 0:
                self.meta["is_temp"] = is_temp
            else:
                raise ValueError("is_temp below absolute zero, not changed")

        self.compute_lamp_emission()
        self.compute_fp_emission()

    def set_bb_aperture(self, value: float):
        """Change the flux-controlling aperture for the black-body source

        Currently, this accepts a float value between 0 and 1 that gives the
        fraction of black-body emission that exits the integrating sphere
        compared to the fully open case.
        """
        if value > 1 or value < 0:
            value = max(0, min(value, 1))
            logger.warning("bb_aperture value out of range [0, 1], clipping to %f",
                           value)
        self.bb_aperture = value
        self.compute_lamp_emission()

    def set_fpmask(self, maskname: str = None, angle: float = None, shift: tuple = None):
        """Change the focal-plane mask

        If `maskname` is not given, the currently inserted mask is rotated to `angle`
        and shifted to `shift`. If `maskname` is given, angle and shift are reset to
        zero, unless explicitely specified.

        See also :class:`FPMask`.

        Parameters
        ----------
        maskname: str, Path
            Name of the mask, as a filepath or to be resolved in irdb
        angle: float [deg]
            Angle by which mask is rotated
        shift: tuple (float, float) [arcsec]
            Shift of mask in x and y direction
        """
        if maskname is not None:
            # Mask is changed: Reset angle and shift
            self.meta["current_fpmask"] = maskname
            self.meta["fpmask_angle"] = 0
            self.meta["fpmask_shift"] = (0, 0)
        if angle is not None:
            self.meta["fpmask_angle"] = angle
        if shift is not None:
            self.meta["fpmask_shift"] = shift

        self.fpmask = FPMask(maskname=self.meta["current_fpmask"],
                             fpmask_filename_format=self.meta['fpmask_filename_format'],
                             angle=self.meta["fpmask_angle"],
                             shift=self.meta["fpmask_shift"])


    @property
    def current_fpmask(self):
        """Name of the mask currently in use"""
        return self.meta['current_fpmask']


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
        with fits.open(path) as hdul:
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
            self.flux_lamp = (self.emiss_bb * self.is_lamp(lam) * self.bb_aperture
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
        """Compute the intensity for the single-line lasers

        The function computes both lasers at once. This is possible because the
        lines are so far apart that there is no band that sees them both.
        """
        lam = self.wavelength
        dlam = lam[1] - lam[0]

        # Laser 1 (L band)     ## TODO move to yaml
        lamc_l = 3.39 * u.um
        power_l = 5e-3 * u.W / (c.c * c.h / lamc_l) * u.ph
        # Laser 2 (tunable), power divided among multiple lines
        lam_t = seq(4.68, 4.78, 0.01) * u.um
        nline = len(lam_t)
        power_t = 70e-3 * u.W / (c.c * c.h / lam) * u.ph / nline
        # Laser 3 (M band)
        lamc_m = 5.26 * u.um
        power_m = 20e-3 * u.W / (c.c * c.h / lamc_m) * u.ph

        sigma = 2 * dlam
        amp = 1/(sigma * np.sqrt(2 * np.pi))

        line_l = Gaussian1D(amplitude=amp, mean=lamc_l, stddev=sigma)
        line_m = Gaussian1D(amplitude=amp, mean=lamc_m, stddev=sigma)
        list_t = [Gaussian1D(amplitude=amp, mean=ll, stddev=sigma) for ll in lam_t]
        line_t = sum(list_t[1:], start=list_t[0])
        flux = (power_l * line_l(lam) + power_m * line_m(lam) + power_t * line_t(lam)) \
            / (np.pi * self.d_is_in**2 / 4)
        intens = flux / (np.pi * u.sr)
        return intens.to(self.bb_scale)

    def compute_fp_emission(self):
        """Compute the emission spectrum from the opaque part of the focal-plane mask"""
        self.wcu_temp = self.meta["wcu_temp"] << u.K

        lam = self.wavelength

        # We assume that T_mask = T_WCU, so that we use RvB's eq.(17) rather than (18)
        # This is independent of the mask emissivity and gives maximum mask emission.
        self.mask_em = BlackBody(self.wcu_temp, scale=self.bb_scale)
        self.intens_fp = self.mask_em(lam)

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
    BlackBody temperature:   {self.meta['bb_temp'] << u.K}
    Integrating sphere temp: {self.meta['is_temp'] << u.K}
    WCU temperature:         {self.meta['wcu_temp'] << u.K}
    BlackBody aperture:      {self.bb_aperture}
    Focal-plane mask:        {self.meta['current_fpmask']}  {self.meta['fpmasks']}"""



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
