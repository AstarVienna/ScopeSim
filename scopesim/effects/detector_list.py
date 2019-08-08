import numpy as np
from astropy import units as u

from .effects import Effect
from .apertures import ApertureMask
from .. import utils
from ..optics.image_plane_utils import header_from_list_of_xy, calc_footprint

__all__ = ["DetectorList"]


class DetectorList(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = [90, 120, 500]
        self.meta["pixel_scale"] = "!INST.pixel_scale"      # arcsec
        self.meta.update(kwargs)

    def fov_grid(self, which="edges", **kwargs):
        """Returns an ApertureMask object. kwargs are "pixel_scale" [arcsec]"""
        aperture_mask = None
        if which == "edges":
            self.meta.update(kwargs)
            self.meta = utils.from_currsys(self.meta)

            hdr = self.image_plane_header
            x_mm, y_mm = calc_footprint(hdr, "D")
            pixel_size = hdr["CDELT1D"]              # mm
            pixel_scale = self.meta["pixel_scale"]   # ["]
            x_sky = x_mm * pixel_scale / pixel_size  # x["] = x[mm] * ["] / [mm]
            y_sky = y_mm * pixel_scale / pixel_size  # y["] = y[mm] * ["] / [mm]

            aperture_mask = ApertureMask(array_dict={"x": x_sky, "y": y_sky},
                                         pixel_scale=pixel_scale)

        return aperture_mask

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

        if ids is None:
            ids = range(len(self.table))

        hdrs = []
        for ii in ids:
            row = self.table[ii]
            xcen, ycen = row["x_cen"], row["y_cen"]
            dx, dy = row["xhw"], row["yhw"]
            cdelt = row["pixsize"]

            hdr = header_from_list_of_xy([xcen-dx, xcen+dx], [ycen-dy, ycen+dy],
                                         pixel_scale=cdelt, wcs_suffix="D")
            if abs(row["angle"]) > 1E-4:
                sang = np.sin(row["angle"] / 57.29578)
                cang = np.cos(row["angle"] / 57.29578)
                hdr["PC1_1"] = cang
                hdr["PC1_2"] = sang
                hdr["PC2_1"] = -sang
                hdr["PC2_2"] = cang

            # hdr["GAIN"] = row["gain"]
            # if "id" in row:
            #     hdr["ID"] = row["id"]
            row_dict = {col: row[col] for col in row.colnames}
            hdr.update(row_dict)
            hdrs += [hdr]

        return hdrs

    def detector_row_dicts(self, ids=None):
        if ids is None:
            ids = range(len(self.table))

        row_dicts = []
        for ii in ids:
            row = self.table[ii]
            row_dicts += [{col: row[col] for col in row.colnames}]

        return row_dicts


def sky_hdr_from_detector_hdr(header, pixel_scale):
    """ pixel_scale in degrees - returns header"""

    x_mm, y_mm = calc_footprint(header, "D")
    scale = pixel_scale / header["CDELT1D"]          # (deg pix-1) / (mm pix-1)
    header = header_from_list_of_xy(x_mm * scale, y_mm * scale, pixel_scale)

    return header
