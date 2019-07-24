import numpy as np

from .effects import Effect
from ..utils import zendist2airmass, airmass2zendist, from_currsys

from .. import rc


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)
        self.meta["z_order"] = [30, 330]

    def apply_to(self, obj, **kwargs):
        return obj

    def fov_grid(self, which="shifts", **kwargs):
        waves, dx, dy = [], [], []
        return [waves, dx, dy]


class AtmosphericDispersion(Shift3D):
    """
    Used to generate the wavelength bins based on shifts due to the atmosphere

    Doesn't contain an ``apply_to`` function, but provides information through
    the ``fov_grid`` function.

    Required Parameters
    -------------------
    airmass : float
        Recommended to use "!OBS.airmass" in the OBS properties
    temperature : float
        [degC] Recommended to use "!ATMO.temperature" in the ATMO properties
    humidity : float
        [0..1] Recommended to use "!ATMO.humidity" in the ATMO properties
    pressure : float
        [bar] Recommended to use "!ATMO.pressure" in the ATMO properties
    latitude : float
        [deg] Recommended to use "!ATMO.latitude" in the ATMO properties
    altitude
        [m] Recommended to use "!ATMO.altitude" in the ATMO properties
    pixel_scale
        [arcsec] Recommended to use "!INST.pixel_scale" in the INST properties

    Optional Parameters
    -------------------
    wave_min : float
        [um] Defaults to "!SIM.spectral.lam_min"
    wave_mid : float
        [um] Defaults to "!SIM.spectral.lam_mid"
    wave_max : float
        [um] Defaults to "!SIM.spectral.lam_max"
    sub_pixel_fraction : float
        [0..1] Defaults to "!SIM.sub_pixel.fraction"
    num_steps : int
        Default: 1000. Number of wavelength steps to use when interpolating the
        atmospheric dispersion curve

    """
    def __init__(self, **kwargs):
        super(AtmosphericDispersion, self).__init__(**kwargs)
        self.meta["z_order"] = [31, 331]
        self.meta["wave_min"] = "!SIM.spectral.lam_min"
        self.meta["wave_mid"] = "!SIM.spectral.lam_mid"
        self.meta["wave_max"] = "!SIM.spectral.lam_max"
        self.meta["sub_pixel_fraction"] = "!SIM.sub_pixel.fraction"
        self.meta["num_steps"] = 1000

        required_keys = ["airmass", "temperature", "humidity", "pressure",
                         "latitude", "altitude", "pupil_angle", "pixel_scale"]
        if not all([key in self.meta for key in required_keys]):
            raise ValueError("One or more of the following keys missing from "
                             "self.meta: \n{} \n{}"
                             "".format(required_keys, self.meta.keys()))

    def fov_grid(self, which="shifts", **kwargs):
        """
        Notes
        -----
        Success! Returns the same values as:
        http://gtc-phase2.gtc.iac.es/science/astroweb/atmosRefraction.php

        """

        for key in self.meta:
            self.meta[key] = from_currsys(self.meta[key])

        atmo_params = {"z0"     : airmass2zendist(self.meta["airmass"]),
                       "temp"   : self.meta["temperature"],         # in degC
                       "rel_hum": self.meta["humidity"] * 100,      # in %
                       "pres"   : self.meta["pressure"] * 1000,     # in mbar
                       "lat"    : self.meta["latitude"],            # in deg
                       "h"      : self.meta["altitude"]}            # in m
        self.meta.update(atmo_params)

        waves, shifts = get_pixel_border_waves_from_atmo_disp(**self.meta)
        dx = shifts * np.cos(np.deg2rad(self.meta["pupil_angle"]))
        dy = shifts * np.sin(np.deg2rad(self.meta["pupil_angle"]))

        if which == "shifts":
            return [waves, dx, dy]
        else:
            return None


class AtmosphericDispersionCorrection(Shift3D):
    def __init__(self, **kwargs):
        """
        Alters the position on the detector for a FOV object (WCS_prefix="D")
        Only acts on FOVs during the main effects loop in OpticalTrain

        Parameters
        ----------
        kwargs
        """
        super(AtmosphericDispersionCorrection, self).__init__(**kwargs)
        self.meta["z_order"] = [2, 302]

    def apply_to(self, obj, **kwargs):
        self.meta.update(kwargs)
        airmass = self.meta["airmass"]
        efficiency = self.meta["efficiency"] if "efficiency" in self.meta else 1

        if airmass == "OBS_AIRMASS":
            # use the same as the atmosphere uses
            pass

        return obj


def atmospheric_refraction(lam, z0=60, temp=0, rel_hum=60, pres=750,
                           lat=-24.5, h=3064):
    """Compute atmospheric refraction

    The function computes the angular difference between the apparent position
    of a star seen from the ground and its true position.

    Parameters
    ----------
    lam : float, np.ndarray
        [um] wavelength bin centres
    z0 : float, optional
        [deg] zenith distance. Default is 60 deg from zenith
    temp : float, optional
        [deg C] ground temperature. Default is 0 deg C
    rel_hum : float, optional
        [%] relative humidity. Default is 60%
    pres : float, optional
        [millibar] air pressure. Default is 750 mbar
    lat : float, optional
        [deg] latitude. Default set for Cerro Armazones: 24.5 deg South
    h : float, optional
        [m] height above sea level. Default is 3064 m

    Returns
    -------
    ang : float, np.ndarray
        [arcsec] angle between real position and refracted position

    References
    ----------
    See Stone 1996 and the review by S. Pedraz -
    http://www.caha.es/newsletter/news03b/pedraz/newslet.html
    """

    # need T, P, RH for Ps, Pw Pa
    T = 273.15 + temp

    Ps = -10474. + 116.43 * T - 0.43284 * T**2 + 0.0005384 * T**3
    Pw = Ps * rel_hum / 100.
    Pa = pres - Pw

    # need n0 for gamma
    sig = 1. / lam
    Da = (Pa / T) * (1. + Pa * (57.9E-8 - 0.0009325 / T + 0.25844 / T**2))
    Dw = (Pw / T) * (1. + Pw * (1. + 3.7E-4 * Pw) *
                     (-2.37321E-3 + 2.23366 / T - 710.792 / T**2 +
                      77514.1 / T**3))

    na = Da * (2371.34 + 683939.7 / (130. - sig**2) + 4547.3 / (38.9 - sig**2))
    nw = Dw * (6487.31 + 58.058 * sig**2 - 0.7115 * sig**4 + 0.08851 * sig**6)
    n0 = 1E-8 * (na + nw) + 1.

    # need gamma, kappa and beta for R
    g = n0 - 1.
    b = 0.001254 * (273.15 + temp) / 273.15
    k = 1. + 0.005302 * np.sin(np.deg2rad(lat))**2 \
        - 0.00000583 * np.sin(2. * np.deg2rad(lat))**2 - 0.000000315 * h

    R = k * g * (1 - b) * np.tan(np.deg2rad(z0)) \
        - k * g * (b - g/2.) * np.tan(np.deg2rad(z0))**3

    # the refraction is the side of a triangle, although the triangle
    # side makes an arc across the sky.
    # We want the angle that this triangle side is subtending
    # Using the small angle approximation (which is in radians),
    # we can get the angle of refraction

    ang = np.rad2deg(R) * 3600

    # return value is in arcsec
    return ang


def get_pixel_border_waves_from_atmo_disp(**kwargs):
    atmo_disp_dict = {key: kwargs[key] for key in ["z0", "temp", "rel_hum",
                                                   "pres", "lat", "h"]}
    wave_range = np.linspace(kwargs["wave_min"], kwargs["wave_max"],
                             kwargs["num_steps"])
    wave_mid = kwargs["wave_mid"]
    offset_mid = atmospheric_refraction(lam=wave_mid, **atmo_disp_dict)
    offset_ang = atmospheric_refraction(lam=wave_range, **atmo_disp_dict)
    offset_ang -= offset_mid

    offset_step = kwargs["pixel_scale"] * kwargs["sub_pixel_fraction"]
    offset_pix = offset_ang / offset_step

    # interpolate to get the edge wavelengths of the pixels
    y = wave_range[::-1]
    x = offset_pix[::-1]
    xnew = np.unique(x.astype(int))
    xnew = np.array([xnew[0] - 1] + list(xnew) + [xnew[-1] + 1])
    ynew = np.interp(xnew, x, y)

    shifts_angle_edges = xnew[::-1] * offset_step
    wave_pixel_edges = ynew[::-1]

    return wave_pixel_edges, shifts_angle_edges
