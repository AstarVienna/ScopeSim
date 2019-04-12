from astropy.io import fits

from ..image_plane_utils import calc_footprint
from .effects import Effect


class ApertureList(Effect):
    def __init__(self, **kwargs):
        super(ApertureList, self).__init__(**kwargs)


class ApertureMask(Effect):
    def __init__(self, **kwargs):
        super(ApertureMask, self).__init__(**kwargs)
        self.meta["z_order"] = [0]

    def fov_grid(self, header=None, waverange=None, **kwargs):
        x_sky, y_sky = calc_footprint(self.header)
        return {"wavelengths": None, "edges": [x_sky, y_sky]}

    @property
    def header(self):
        return fits.Header()


class SpectralTraceList(Effect):
    def __init__(self, **kwargs):
        super(SpectralTraceList, self).__init__(**kwargs)
        self.meta["z_order"] = [0]
