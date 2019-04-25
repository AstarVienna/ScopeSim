from scopesim.optics.surface import SpectralSurface
from .effects import Effect


class TERCurve(Effect):
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)

        self.surface = SpectralSurface()
        self.surface.meta.update(self.meta)
        data = self.get_data()
        if data is not None:
            self.surface.table = data


class AtmosphericTERCurve(TERCurve):
    def __init__(self, **kwargs):
        super(AtmosphericTERCurve, self).__init__(**kwargs)


class SkycalcTERCurve(TERCurve):
    def __init__(self, **kwargs):
        super(SkycalcTERCurve, self).__init__(**kwargs)


class QuantumEfficiencyCurve(TERCurve):
    def __init__(self, **kwargs):
        super(QuantumEfficiencyCurve, self).__init__(**kwargs)
