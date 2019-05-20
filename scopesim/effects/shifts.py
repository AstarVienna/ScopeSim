from .effects import Effect


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)
        self.meta["z_order"] = [0]

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
        """


class AtmosphericDispersionCorrection(Shift3D):
    def __init__(self, **kwargs):
        super(AtmosphericDispersionCorrection, self).__init__(**kwargs)

    def apply_to(self, obj, **kwargs):
        self.meta.update(kwargs)
        airmass = self.meta["airmass"]
        efficiency = self.meta["efficiency"] if "efficiency" in self.meta else 1

        if airmass == "OBS_AIRMASS":
            # use the same as the atmosphere uses
            pass
