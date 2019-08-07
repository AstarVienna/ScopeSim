from .effects import Effect


class SpectralTraceList(Effect):
    def __init__(self, **kwargs):
        super(SpectralTraceList, self).__init__(**kwargs)
        self.meta["z_order"] = [70]
