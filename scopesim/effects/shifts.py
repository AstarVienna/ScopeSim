import numpy as np

from .effects import Effect
from .effects_utils import atmospheric_refraction as atmo_refr
from ..utils import zendist2airmass, airmass2zendist

from .. import rc


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]

    def apply_to(self, obj, **kwargs):
        return obj

    def fov_grid(self, which="shifts", **kwargs):
        waves, dx, dy = [], [], []
        return [waves, dx, dy]


class AtmosphericDispersion(Shift3D):
    def __init__(self, **kwargs):
        super(AtmosphericDispersion, self).__init__(**kwargs)

        """
        Needed parameters from atmospheric optical element
        altitude
        latitude
        airmass
        temperature
        humidity
        pressure
        fov_grid
        
        Alters the position on the sky for a FOV object (WCS_prefix="")
        Only acts on FOVs when FOVs are initialised in FOV_Manager
        
        """
        self.meta["z_order"] = [1, 301]

        required_keys = ["airmass", "temperature", "humidity", "pressure",
                         "latitude", "altitude", "pupil_angle"]
        if not all([key in self.meta for key in required_keys]):
            raise ValueError("One or more of the following keys missing from "
                             "self.meta: \n{} \n{}"
                             "".format(required_keys, self.meta.keys()))

    def fov_grid(self, which="shifts", **kwargs):
        """

        Parameters
        ----------
        which
        kwargs

        Returns
        -------

        Notes
        -----
        Success! Returns the same values as:
        http://gtc-phase2.gtc.iac.es/science/astroweb/atmosRefraction.php

        """

        lam_min = rc.__old_config__["SIM_LAM_MIN"]
        lam_mid = rc.__old_config__["SIM_LAM_MID"]
        lam_max = rc.__old_config__["SIM_LAM_MAX"]

        atmo_params = {"z0"     : airmass2zendist(self.meta["airmass"]),
                       "temp"   : self.meta["temperature"],         # in degC
                       "rel_hum": self.meta["humidity"] * 100,      # in %
                       "pres"   : self.meta["pressure"] * 1000,     # in mbar
                       "lat"    : self.meta["latitude"],
                       "h"      : self.meta["altitude"]}

        offset_mid = atmo_refr(lam_mid, **atmo_params)
        offset_min, offset_max = atmo_refr(np.array([lam_min, lam_max]),
                                           **atmo_params) - offset_mid



        waves, dx, dy = [], [], []
        return [waves, dx, dy]


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
