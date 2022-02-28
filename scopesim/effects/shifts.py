import numpy as np
from astropy import units as u
from astropy.table import Table

from .effects import Effect
from ..utils import airmass2zendist, from_currsys, check_keys, quantify
from ..base_classes import FieldOfViewBase


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)
        params = {"z_order": [30, 230],
                  "report_plot_include": True,
                  "report_table_include": False,}
        self.meta.update(params)
        self.meta.update(kwargs)

    def apply_to(self, obj, **kwargs):
        return obj

    def fov_grid(self, which="shifts", **kwargs):
        if which == "shifts":
            col_names = ["wavelength", "dx", "dy"]
            waves, dx, dy = [self.get_table(**kwargs)[col] for col in col_names]
            return waves, dx, dy
        else:
            return None

    def get_table(self, **kwargs):
        if self.table is None:
            names = ["wavelength", "dx", "dy"]
            waves = from_currsys(["!SIM.spectral.wave_" + key
                                  for key in ("min", "mid", "max")])
            tbl = Table(names=names, data=[waves, [0] * 3, [0] * 3])
        else:
            tbl = self.table

        return tbl

    def plot(self):
        import matplotlib.pyplot as plt
        plt.gcf().clf()

        tbl = self.get_table()
        plt.scatter(x=tbl["dx"], y=tbl["dy"], c=tbl["wavelength"])
        plt.colorbar()
        plt.xlabel("dx [{}]".format(quantify(tbl["dx"], u.arcsec).unit))
        plt.ylabel("dy [{}]".format(quantify(tbl["dy"], u.arcsec).unit))
        plt.axvline(0, ls=":")
        plt.axhline(0, ls=":")
        # plt.gca().set_aspect("equal")

        return plt.gcf()


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
        [um] Defaults to "!SIM.spectral.wave_min"
    wave_mid : float
        [um] Defaults to "!SIM.spectral.wave_mid"
    wave_max : float
        [um] Defaults to "!SIM.spectral.wave_max"
    sub_pixel_fraction : float
        [0..1] Defaults to "!SIM.sub_pixel.fraction"
    num_steps : int
        Default: 1000. Number of wavelength steps to use when interpolating the
        atmospheric dispersion curve

    """
    def __init__(self, **kwargs):
        super(AtmosphericDispersion, self).__init__(**kwargs)
        params = {"z_order": [231],
                  "wave_min": "!SIM.spectral.wave_min",
                  "wave_mid": "!SIM.spectral.wave_mid",
                  "wave_max": "!SIM.spectral.wave_max",
                  "sub_pixel_fraction": "!SIM.sub_pixel.fraction",
                  "num_steps": 1000,}
        self.meta.update(params)
        self.meta.update(kwargs)

        required_keys = ["airmass", "temperature", "humidity", "pressure",
                         "latitude", "altitude", "pupil_angle", "pixel_scale"]
        check_keys(self.meta, required_keys, action="error")

    def get_table(self, **kwargs):
        """
        Called by the fov_grid method of Shift3D

        Returns
        -------
        tbl : astropy.Table
            A table with the columns
            - waves : [um]
            - dx, dy : [arcsec]

        Notes
        -----
        Success! Returns the same values as:
        http://gtc-phase2.gtc.iac.es/science/astroweb/atmosRefraction.php
        """

        if len(kwargs) > 0:
            self.meta.update(kwargs)

        airmass = from_currsys(self.meta["airmass"])
        atmo_params = {"z0": airmass2zendist(airmass),
                       "temp": self.meta["temperature"],  # in degC
                       "rel_hum": self.meta["humidity"] * 100,  # in %
                       "pres": self.meta["pressure"] * 1000,  # in mbar
                       "lat": self.meta["latitude"],  # in deg
                       "h": self.meta["altitude"]}  # in m
        self.meta.update(atmo_params)
        params = from_currsys(self.meta)

        waves, shifts = get_pixel_border_waves_from_atmo_disp(**params)
        dx = shifts * np.sin(np.deg2rad(params["pupil_angle"]))
        dy = shifts * np.cos(np.deg2rad(params["pupil_angle"]))

        names = ["wavelength", "dx", "dy"]
        tbl = Table(names=names, data=[waves, dx, dy])

        return tbl


class AtmosphericDispersionCorrection(Shift3D):
    def __init__(self, **kwargs):
        """
        Alters the position on the detector for a FOV object (WCS_prefix="D")

        Only acts on FOVs during the main effects loop in OpticalTrain.
        For the sake of computational efficiency, the ADC can be instructed to
        counteract the atmospheric diffraction during the OpticalTrain setup
        phase, by passing the kwarg: ``quick_adc=True``

        Parameters
        ----------
        kwargs
        """
        super(AtmosphericDispersionCorrection, self).__init__(**kwargs)
        self.meta["z_order"] = [632]
        if "quick_adc" in self.meta and self.meta["quick_adc"] is True:
            self.meta["z_order"] += [232]
        if "efficiency" not in self.meta:
            self.meta["efficiency"] = 1
        self.apply_to_classes = FieldOfViewBase

        required_keys = ["airmass", "temperature", "humidity", "pressure",
                         "latitude", "altitude", "pupil_angle", "pixel_scale",
                         "wave_mid"]
        check_keys(self.meta, required_keys, action="error")

        if self.table is None:
            self.table = self.get_table()

    def apply_to(self, fov, **kwargs):
        # .. todo:: Currently applying shift with pixel_scale to CRPIX-D
        # .. todo:: Change this to be applying to CRVAL-D using plate_scale
        # get mid wavelength of fov
        # work out on-sky shift
        # work out pupil-plane shift
        # convert to pixel shifts
        # correct fovs CRPIXnD keys

        if isinstance(fov, self.apply_to_classes):
            self.meta = from_currsys(self.meta)
            atmo_params = {"z0":        airmass2zendist(self.meta["airmass"]),
                           "temp":      self.meta["temperature"],  # in degC
                           "rel_hum":   self.meta["humidity"] * 100,  # in %
                           "pres":      self.meta["pressure"] * 1000, # in mbar
                           "lat":       self.meta["latitude"],  # in deg
                           "h":         self.meta["altitude"]}  # in m
            fov_wave_mid = quantify(fov.wavelength, u.um).value
            # .. todo:: why are we comparing to shift_obs. Is this static?
            obs_wave_mid = self.meta["wave_mid"]
            shift_obs = atmospheric_refraction(lam=obs_wave_mid, **atmo_params)
            shift_fov = atmospheric_refraction(lam=fov_wave_mid, **atmo_params)
            shift_rel_arcsec = shift_fov - shift_obs
            shift_rel_pixel = shift_rel_arcsec / self.meta["pixel_scale"]

            # this assumes a one-to-one mapping of sky coords to detector coords
            # this won't work if there is appreciable distortion in the optics
            # .. todo:: think about this
            pup_ang = self.meta["pupil_angle"]
            dy_pix = shift_rel_pixel * np.cos(np.deg2rad(pup_ang))
            dx_pix = shift_rel_pixel * np.sin(np.deg2rad(pup_ang))

            fov.header["CRPIX1D"] += dx_pix
            fov.header["CRPIX2D"] += dy_pix

        return fov

    def fov_grid(self, which="shifts", **kwargs):
        kwargs.update(self.meta)
        if "quick_adc" in self.meta:
            ad = AtmosphericDispersion(**self.meta)
            waves, dx, dy = ad.fov_grid()
            dx *= -(1 - self.meta["efficiency"])
            dy *= -(1 - self.meta["efficiency"])
            return waves, dx, dy
        else:
            return None

    def plot(self):
        return None


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
    """

    Parameters
    ----------
    kwargs

    Returns
    -------

    """
    atmo_disp_dict = {key: kwargs[key] for key in ["z0", "temp", "rel_hum",
                                                   "pres", "lat", "h"]}
    wave_range = np.linspace(kwargs["wave_min"], kwargs["wave_max"],
                             kwargs["num_steps"])
    wave_mid = kwargs["wave_mid"]
    offset_mid = atmospheric_refraction(lam=wave_mid, **atmo_disp_dict)
    offset_ang = atmospheric_refraction(lam=wave_range, **atmo_disp_dict)
    offset_ang -= offset_mid

    # .. todo:: replace the 1e-7 with a variable in !SIM
    if np.any(np.abs(offset_ang) > 1e-7):
        offset_step = kwargs["pixel_scale"] * kwargs["sub_pixel_fraction"]
        offset_pix = offset_ang / offset_step

        # interpolate to get the edge wavelengths of the pixels
        # off_new is always increasing, thanks to np.unique
        off_new = np.unique(offset_pix.astype(int))
        # add 1 pixel to either end of the range covered by (wave_min, wave_max)
        off_new = np.array([off_new[0]-1] + list(off_new) + [off_new[-1]+1])

        if offset_pix[0] > offset_pix[-1]:
            wave_range = wave_range[::-1]
            offset_pix = offset_pix[::-1]

        wave_new = np.interp(off_new, offset_pix, wave_range)

        if wave_new[0] > wave_new[-1]:
            wave_new = wave_new[::-1]
            off_new = off_new[::-1]

        shifts_angle_edges = off_new * offset_step
        wave_pixel_edges = wave_new

    else:
        wave_pixel_edges = np.array([kwargs["wave_min"], kwargs["wave_max"]])
        shifts_angle_edges = np.zeros(2)

    return wave_pixel_edges, shifts_angle_edges
