from .effects import Effect


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 300]

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
        """
        self.meta["z_order"] = [1, 301]

    def apply_to(self, obj, **kwargs):
        return obj

    def fov_grid(self, which="shifts", **kwargs):
        waves, dx, dy = [], [], []
        return [waves, dx, dy]





class AtmosphericDispersionCorrection(Shift3D):
    def __init__(self, **kwargs):
        """
        Alters the position on the detector for a FOV object (WCS_prefix="D")

        Parameters
        ----------
        kwargs
        """
        super(AtmosphericDispersionCorrection, self).__init__(**kwargs)

    def apply_to(self, obj, **kwargs):
        self.meta.update(kwargs)
        airmass = self.meta["airmass"]
        efficiency = self.meta["efficiency"] if "efficiency" in self.meta else 1

        if airmass == "OBS_AIRMASS":
            # use the same as the atmosphere uses
            pass

        return obj
