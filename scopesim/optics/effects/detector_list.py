import numpy as np
from astropy import units as u

from ... import utils
from .effects import Effect
from ..image_plane_utils import header_from_list_of_xy, calc_footprint


class DetectorList(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [0, 100, 500]

    def fov_grid(self, header=None, waverange=None, **kwargs):
        hdr = self.image_plane_header
        pixel_scale = kwargs["pixel_scale"]  # ["] --> [deg]
        pixel_size = hdr["CDELT1D"]
        x_mm, y_mm = calc_footprint(hdr, "D")
        x_sky = x_mm * pixel_scale / pixel_size  # x[deg] = x[mm] * [deg] / [mm]
        y_sky = y_mm * pixel_scale / pixel_size  # y[deg] = y[mm] * [deg] / [mm]

        return {"edges": [x_sky, y_sky], "wavelengths": None}

    @property
    def image_plane_header(self):
        tbl = self.get_data()
        pixel_size = np.min(utils.quantity_from_table("pixsize", tbl, u.mm))
        x_unit = utils.unit_from_table("x_cen", tbl, u.mm)
        y_unit = utils.unit_from_table("y_cen", tbl, u.mm)

        x_det_min = np.min(tbl["x_cen"] - tbl["xhw"]) * x_unit
        x_det_max = np.max(tbl["x_cen"] + tbl["xhw"]) * x_unit
        y_det_min = np.min(tbl["y_cen"] - tbl["yhw"]) * y_unit
        y_det_max = np.max(tbl["y_cen"] + tbl["yhw"]) * y_unit

        x_det = [x_det_min.to(u.mm).value, x_det_max.to(u.mm).value]
        y_det = [y_det_min.to(u.mm).value, y_det_max.to(u.mm).value]

        pixel_size = pixel_size.to(u.mm).value
        hdr = header_from_list_of_xy(x_det, y_det, pixel_size, "D")

        return hdr

    def detector_headers(self, ids=None):
        """Return the header to describe each individual detector"""
        raise NotImplementedError


def sky_hdr_from_detector_hdr(header, pixel_scale):
    """ pixel_scale in degrees - returns header"""

    pixel_size = header["CDELT1D"]
    x_mm, y_mm = calc_footprint(header, "D")
    scale = pixel_scale / pixel_size            # (deg pix-1) / (mm pix-1)
    header = header_from_list_of_xy(x_mm * scale, y_mm * scale, pixel_scale)

    return header
