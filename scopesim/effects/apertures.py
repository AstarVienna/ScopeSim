import numpy as np
from matplotlib.path import Path
from astropy.io import fits
from astropy import units as u

from .effects import Effect
from ..optics import image_plane_utils as imp_utils

from ..utils import quantity_from_table, from_currsys


class ApertureMask(Effect):
    """
    Only provides the on-sky window coords of the Aperture

    - Case: Imaging
        - Covers the whole FOV of the detector
        - Round (with mask), square (without mask)
    - Case : LS Spec
        - Covers the slit FOV
        - Polygonal (with mask), square (without mask)
    - Case : IFU Spec
        - Covers the on-sky FOV of one slice of the IFU
        - Square (without mask)
    - Case : MOS Spec
        - Covers a single MOS fibre FOV
        - Round, Polygonal (with mask), square (without mask)

    """
    def __init__(self, **kwargs):
        super(ApertureMask, self).__init__(**kwargs)
        self.meta["z_order"] = [80, 110]
        self.meta["pixel_scale"] = "!INST.pixel_scale"
        self.meta["no_mask"] = True
        self.meta["angle"] = 0

        self.meta.update(kwargs)

        self._header = None
        self._mask = None
        self.mask_sum = None

    def fov_grid(self, which="edges", **kwargs):
        """ Returns a header with the sky coordinates """
        self.meta.update(kwargs)
        return self.header

    @property
    def hdu(self):
        return fits.ImageHDU(data=self.mask, header=self.header)

    @property
    def header(self):
        if not isinstance(self._header, fits.Header) \
                and "x" in self.table.colnames and "y" in self.table.colnames:
            self._header = self.get_header()
        return self._header

    def get_header(self):
        self.meta = from_currsys(self.meta)
        x = quantity_from_table("x", self.table, u.arcsec).to(u.deg).value
        y = quantity_from_table("y", self.table, u.arcsec).to(u.deg).value
        pix_scale_deg = self.meta["pixel_scale"] / 3600.
        header = imp_utils.header_from_list_of_xy(x, y, pix_scale_deg)

        return header

    @property
    def mask(self):
        if not isinstance(self._header, fits.Header) \
                and "x" in self.table.colnames and "y" in self.table.colnames:
            self._mask = self.get_mask()
        return self._mask

    def get_mask(self):
        if self.meta["no_mask"] is False:
            x = quantity_from_table("x", self.table, u.arcsec).to(u.deg).value
            y = quantity_from_table("y", self.table, u.arcsec).to(u.deg).value
            pixel_scale_deg = self.meta["pixel_scale"] / 3600.
            angle_deg = self.meta["angle"]
            mask = mask_from_coords(x, y, angle_deg, pixel_scale_deg)
        else:
            mask = None

        return mask


def mask_from_coords(x0, y0, angle, pixel_scale):
    angle = np.deg2rad(angle)
    c, s = np.cos(angle), np.sin(angle)
    x = c * x0 - s * y0
    y = s * x0 + c * y0

    hdr = imp_utils.header_from_list_of_xy(x, y, pixel_scale)
    ra, dec = imp_utils.calc_footprint(hdr)
    naxis1, naxis2 = hdr["NAXIS1"], hdr["NAXIS2"]
    xrange = np.linspace(np.min(ra), np.max(ra), naxis1)
    yrange = np.linspace(np.min(dec), np.max(dec), naxis2)
    coords = [(xi, yi) for xi in xrange for yi in yrange]

    path = Path([(xi, yi) for xi, yi in zip(x, y)])

    # ..todo:: replace the 0.005 with a variable in !SIM.computing
    mask = path.contains_points(coords, radius=0.005).reshape((naxis1, naxis2))

    return mask


class ApertureList(Effect):
    def __init__(self, **kwargs):
        super(ApertureList, self).__init__(**kwargs)
        self.meta["z_order"] = [81, 111]

    def fov_grid(self, which="edges", **kwargs):
        """ Returns a list of ApertureMasks """
        return []


class SquareApertureList(ApertureList):
    def __init__(self, **kwargs):
        super(SquareApertureList, self).__init__(**kwargs)
        self.meta["z_order"] = [82, 112]


class RoundApertureList(ApertureList):
    def __init__(self, **kwargs):
        super(RoundApertureList, self).__init__(**kwargs)
        self.meta["z_order"] = [83, 113]


class PolygonApertureList(ApertureList):
    def __init__(self, **kwargs):
        super(PolygonApertureList, self).__init__(**kwargs)
        self.meta["z_order"] = [84, 114]
