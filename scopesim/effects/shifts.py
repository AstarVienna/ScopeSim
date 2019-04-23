from .effects import Effect


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)
        self.meta["z_order"] = [0]

    def fov_grid(self, header=None, waverange=None, **kwargs):
        dic = {"wavelengths": waverange, "x_shifts": [0, 0], "y_shifts": [0, 0]}
        return dic


class AtmosphericDispersion(Shift3D):
    def __init__(self, **kwargs):
        super(AtmosphericDispersion, self).__init__(**kwargs)


class AtmosphericDispersionCorrector(Shift3D):
    def __init__(self, **kwargs):
        super(AtmosphericDispersionCorrector, self).__init__(**kwargs)
