'''Effects related to field masks, including spectroscopic slits'''
from os import path as pth
from copy import deepcopy
import logging

import numpy as np
from matplotlib.path import Path
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from .effects import Effect
from ..optics import image_plane_utils as imp_utils
from ..base_classes import FOVSetupBase

from ..utils import quantify, quantity_from_table, from_currsys, check_keys


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

    The geometry of an ``ApertureMask`` can be initialised with the standard
    DataContainer methods (see Parameters below). Regardless of which method
    is used, the following columns must be present::

        x       y
        arcsec  arcsec
        float   float

    Certain keywords need to also be included in the ascii header::

        # id: <int>
        # conserve_image: <bool>
        # x_unit: <str>
        # y_unit: <str>

    If ``conserve_image`` is ``False``, the flux from all sources in the
    aperture is summed and distributed uniformly over the aperture area.


    Parameters
    ----------
    filename : str
        Path to ASCII file containing the columns listed above

    table : astropy.Table
        An astropy Table containing the columns listed above

    array_dict : dict
        A dictionary containing the columns listed above:
        ``{x: [...], y: [...], id: <int>, conserve_image: <bool>}

    Other Parameters
    ----------------
    pixel_scale : float
        [arcsec] Defaults to ``"!INST.pixel_scale"`` from the config

    id : int
        An integer to identify the ``ApertureMask`` in a list of apertures

    """
    def __init__(self, **kwargs):
        if not np.any([key in kwargs for key in ["filename", "table",
                                                 "array_dict"]]):
            if "width" in kwargs and "height" in kwargs and \
                    "filename_format" in kwargs:
                kwargs = from_currsys(kwargs)
                w, h = kwargs["width"], kwargs["height"]
                kwargs["filename"] = kwargs["filename_format"].format(w, h)

        super(ApertureMask, self).__init__(**kwargs)
        params = {"pixel_scale": "!INST.pixel_scale",
                  "no_mask": True,
                  "angle": 0,
                  "shape": "rect",
                  "conserve_image": True,
                  "id": 0,
                  "report_plot_include": False,
                  "report_table_include": True,
                  "report_table_rounding": 4}

        self.meta["z_order"] = [80, 280, 380]
        self.meta.update(params)
        self.meta.update(kwargs)

        self._header = None
        self._mask = None
        self.mask_sum = None

        self.required_keys = ["filename", "table", "array_dict"]
        check_keys(kwargs, self.required_keys, "warning", all_any="any")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, FOVSetupBase):
            x = quantity_from_table("x", self.table, u.arcsec).to(u.arcsec).value
            y = quantity_from_table("y", self.table, u.arcsec).to(u.arcsec).value
            obj.shrink(["x", "y"], ([min(x), max(x)], [min(y), max(y)]))

            # ..todo: HUGE HACK - Get rid of this!
            for vol in obj.volumes:
                vol["meta"]["xi_min"] = min(x) * u.arcsec
                vol["meta"]["xi_max"] = max(x) * u.arcsec

        return obj

    def fov_grid(self, which="edges", **kwargs):
        """ Returns a header with the sky coordinates """
        logging.warning("DetectorList.fov_grid will be depreciated in v1.0")
        if which == "edges":
            self.meta.update(kwargs)
            return self.header
        elif which == "masks":
            self.meta.update(kwargs)
            return self.mask

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
        header["APERTURE"] = self.meta["id"]
        header["ROT"] = self.meta["angle"]
        header["IMG_CONS"] = self.meta["conserve_image"]

        return header

    @property
    def mask(self):
        if not isinstance(self._header, fits.Header) \
                and "x" in self.table.colnames and "y" in self.table.colnames:
            self._mask = self.get_mask()
        return self._mask

    def get_mask(self):
        """
        For placing over FOVs if the Aperture is rotated w.r.t. the field
        """
        self.meta = from_currsys(self.meta)

        if self.meta["no_mask"] is False:
            x = quantity_from_table("x", self.table, u.arcsec).to(u.deg).value
            y = quantity_from_table("y", self.table, u.arcsec).to(u.deg).value
            pixel_scale_deg = self.meta["pixel_scale"] / 3600.
            mask = mask_from_coords(x, y, pixel_scale_deg)
        else:
            mask = None

        return mask

    def plot(self, new_figure=True):
        import matplotlib.pyplot as plt
        if new_figure:
            plt.gcf().clf()

        x = list(self.table["x"].data)
        y = list(self.table["y"].data)
        plt.plot(x + [x[0]], y + [y[0]])
        plt.gca().set_aspect("equal")

        return plt.gcf()


class RectangularApertureMask(ApertureMask):
    def __init__(self, **kwargs):
        super(RectangularApertureMask, self).__init__(**kwargs)
        params = {"x_unit" : "arcsec",
                  "y_unit" : "arcsec"}
        self.meta.update(params)
        self.meta.update(kwargs)
        check_keys(self.meta, ["x", "y", "width", "height"])

        self.table = self.get_table(**kwargs)

    def fov_grid(self, which="edges", **kwargs):
        """ Returns a header with the sky coordinates """
        if which == "edges":
            self.table = self.get_table(**kwargs)
            return self.header      # from base class ApertureMask

        elif which == "masks":
            self.meta.update(kwargs)
            return self.mask

    def get_table(self, **kwargs):
        self.meta.update(kwargs)
        x = from_currsys(self.meta["x"])
        y = from_currsys(self.meta["y"])
        dx = 0.5 * from_currsys(self.meta["width"])
        dy = 0.5 * from_currsys(self.meta["height"])
        xs = [x - dx, x + dx, x + dx, x - dx]
        ys = [y - dy, y - dy, y + dy, y + dy]
        tbl = Table(names=["x", "y"], data=[xs, ys], meta=self.meta)

        return tbl


class ApertureList(Effect):
    """
    A list of apertures, useful for IFU or MOS instruments

    Much like an ApertureMask, an ApertureList can be initialised by either
    of the three standard DataContainer methods. The easiest is however to
    make an ASCII file with the following columns::

        id   left    right   top     bottom  angle  conserve_image  shape
             arcsec  arcsec  arcsec  arcsec  deg
        int  float   float   float   float   float  bool            str/int

    Acceptable ``shape`` entries are: ``round``, ``rect``, ``hex``, ``oct``, or
    an integer describing the number of corners the polygon should have.

    A polygonal mask is generated for a given ``shape`` and will be scaled
    to fit inside the edges of each aperture list row. The corners of each
    aperture defined by shape are found by finding equidistant positions around
    an ellipse constrained by the edges (``left``, ..., etc). An additional
    optional column ``offset`` may be added. This column describes the offset
    from 0 deg to the angle where the first corner is set.

    Additionally, the filename of an ``ApertureMask`` polygon file can be given.
    The geometry of the polygon defined in the file will be scaled to fit
    inside the edges of the row.

    .. note:: ``shape`` values ``"rect"`` and ``4`` do not produce equal results

       Both use 4 equidistant points around an ellipse constrained by
       [``left``, ..., etc]. However ``"rect"`` aims to fill the contraining
       area, while ``4`` simply uses 4 points on the ellipse.
       Consequently, ``4`` results in a diamond shaped mask covering only
       half of the constraining area filled by ``"rect"``.

    """
    def __init__(self, **kwargs):
        super(ApertureList, self).__init__(**kwargs)
        params = {"pixel_scale": "!INST.pixel_scale",
                  "n_round_corners": 32,        # number of corners use to estimate ellipse
                  "no_mask": False,             # .. todo:: is this necessary when we have conserve_image?
                  "report_plot_include": True,
                  "report_table_include": True,
                  "report_table_rounding": 4}
        self.meta["z_order"] = [81, 281]
        self.meta.update(params)
        self.meta.update(kwargs)

        if self.table is not None:
            required_keys = ["id", "left", "right", "top", "bottom", "angle",
                             "conserve_image", "shape"]
            check_keys(self.table.colnames, required_keys)

    def fov_grid(self, which="edges", **kwargs):
        params = deepcopy(self.meta)
        params.update(kwargs)
        if which == "edges":
            return [ap.fov_grid(which=which, **params) for ap in self.apertures]
        if which == "masks":
            return {ap.meta["id"]: ap.mask for ap in self.apertures}

    @property
    def apertures(self):
        return self.get_apertures(range(len(self.table)))

    def get_apertures(self, row_ids):
        if isinstance(row_ids, int):
            row_ids = [row_ids]

        apertures_list = []
        for ii in row_ids:
            row = self.table[ii]
            row_dict = {col: row[col] for col in row.colnames}
            row_dict["n_round"] = self.meta["n_round_corners"]
            array_dict = make_aperture_polygon(**row_dict)
            params = {"id": row["id"],
                      "angle": row["angle"],
                      "shape": row["shape"],
                      "conserve_image": row["conserve_image"],
                      "no_mask": self.meta["no_mask"],
                      "pixel_scale": self.meta["pixel_scale"],
                      "x_unit": "arcsec",
                      "y_unit": "arcsec",
                      "angle_unit": "arcsec"}
            apertures_list += [ApertureMask(array_dict=array_dict, **params)]

        return apertures_list

    def plot(self):
        import matplotlib.pyplot as plt
        plt.gcf().clf()

        for ap in self.apertures:
            ap.plot(new_figure=False)

        return plt.gcf()

    def plot_masks(self):
        import matplotlib.pyplot as plt

        aps = self.apertures
        n = len(aps)
        w = np.ceil(n ** 0.5).astype(int)
        h = np.ceil(n / w).astype(int)

        for ii, ap in enumerate(aps):
            plt.subplot(w, h, ii + 1)
            plt.imshow(ap.mask.T)
        plt.show()

    def __add__(self, other):
        if isinstance(other, ApertureList):
            from astropy.table import vstack
            self.table = vstack([self.table, other.table])

            return self
        else:
            raise ValueError("Secondary argument not of type ApertureList: {}"
                             "".format(type(other)))

    def __getitem__(self, item):
        return self.get_apertures(item)[0]



class SlitWheel(Effect):
    """
    This wheel holds a selection of predefined spectroscopic slits
    and possibly other field masks.

    It should contain an open position.
    A user can define a non-standard slit by directly using the Aperture
    effect.

    Examples
    --------
    ::
        name: slit_wheel
        class: SlitWheel
        kwargs:
            slit_names: []
            filename_format: "MASK_slit_{}.dat
            current_slit: "C"
    """

    def __init__(self, **kwargs):
        required_keys = ["slit_names", "filename_format", "current_slit"]
        check_keys(kwargs, required_keys, action="error")

        super(SlitWheel, self).__init__(**kwargs)

        params = {"z_order": [80, 280, 580],
                  "path": "",
                  "report_plot_include": False,
                  "report_table_include": True,
                  "report_table_rounding": 4}
        self.meta.update(params)
        self.meta.update(kwargs)

        path = pth.join(self.meta["path"],
                        from_currsys(self.meta["filename_format"]))
        self.slits = {}
        for name in self.meta["slit_names"]:
            kwargs["name"] = name
            self.slits[name] = ApertureMask(filename=path.format(name),
                                            **kwargs)

        self.table = self.get_table()


    def apply_to(self, obj, **kwargs):
        '''Use apply_to of current_slit'''
        return self.current_slit.apply_to(obj, **kwargs)


    def fov_grid(self, which="edges", **kwargs):
        return self.current_slit.fov_grid(which=which, **kwargs)

    def change_slit(self, slitname=None):
        '''Change the current slit'''
        if not slitname or slitname in self.slits.keys():
            self.meta['current_slit'] = slitname
            self.include = slitname
        else:
            raise ValueError("Unknown slit requested: " + slitname)

    @property
    def current_slit(self):
        '''Return the currently used slit'''
        currslit = from_currsys(self.meta["current_slit"])
        if not currslit:
            return False
        return self.slits[currslit]

    @property
    def display_name(self):
        return f'{self.meta["name"]} : ' \
               f'[{from_currsys(self.meta["current_slit"])}]'


    def __getattr__(self, item):
        return getattr(self.current_slit, item)

    def get_table(self):
        '''Create a table of slits with centre position, width and length

        Width is defined as the extension in the y-direction, length in the
        x-direction. All values are in milliarcsec.'''
        names = list(self.slits.keys())
        slits = self.slits.values()
        xmax = np.array([slit.data['x'].max() * u.Unit(slit.meta['x_unit'])
                         .to(u.mas) for slit in slits])
        xmin = np.array([slit.data['x'].min() * u.Unit(slit.meta['x_unit'])
                         .to(u.mas) for slit in slits])
        ymax = np.array([slit.data['y'].max() * u.Unit(slit.meta['y_unit'])
                         .to(u.mas) for slit in slits])
        ymin = np.array([slit.data['y'].min() * u.Unit(slit.meta['y_unit'])
                         .to(u.mas) for slit in slits])
        xmax = quantify(xmax, u.mas)
        xmin = quantify(xmin, u.mas)
        ymax = quantify(ymax, u.mas)
        ymin = quantify(ymin, u.mas)

        lengths = xmax - xmin
        widths = ymax - ymin
        x_centres = (xmax + xmin) / 2
        y_centres = (ymax + ymin) / 2
        tbl = Table(names=["name", "x_centre", "y_centre", "length", "width"],
                    data=[names, x_centres, y_centres, lengths, widths])
        return tbl


################################################################################


def make_aperture_polygon(left, right, top, bottom, angle, shape, **kwargs):

    n_round = kwargs["n_round"] if "n_round" in kwargs else 32
    offset = kwargs["offset"] if "offset" in kwargs else 0.

    n_corners = {"rect": 4, "hex": 6, "oct": 8, "round": n_round}
    try:
        shape = int(float(shape))
        n_corners[shape] = shape
    except:
        pass

    x0, y0 = 0.5 * (right + left), 0.5 * (top + bottom)
    dx, dy = 0.5 * (right - left), 0.5 * (top - bottom)
    n = n_corners[shape]

    if isinstance(shape, str) and "rect" in shape:
        dx *= 1.41421356
        dy *= 1.41421356
        offset += 45.

    x, y = points_on_a_circle(n=n, x0=x0, y0=y0, dx=dx, dy=dy, offset=offset)
    if angle != 0.:
        x, y = rotate(x=x, y=y, x0=np.average(x), y0=np.average(y), angle=angle)

    return {"x": x, "y": y}


def points_on_a_circle(n, x0=0, y0=0, dx=1, dy=1, offset=0):
    deg2rad = np.pi / 180
    d_angle = np.arange(0, 360, 360 / n) + offset
    x = x0 + dx * np.cos(d_angle * deg2rad)
    y = y0 + dy * np.sin(d_angle * deg2rad)

    return x, y


def mask_from_coords(x, y, pixel_scale):
    naxis1 = int(np.ceil((np.max(x) - np.min(x)) / pixel_scale))
    naxis2 = int(np.ceil((np.max(y) - np.min(y)) / pixel_scale))
    xrange = np.linspace(np.min(x), np.max(x), naxis1)
    yrange = np.linspace(np.min(y), np.max(y), naxis2)
    coords = [(xi, yi) for xi in xrange for yi in yrange]

    corners = [(xi, yi) for xi, yi in zip(x, y)]
    path = Path(corners)
    # ..todo: known issue - for super thin apertures, the first row is masked
    # rad = 0.005
    rad = 0  # increase this to include slightly more points within the polygon
    mask = path.contains_points(coords, radius=rad).reshape((naxis2, naxis1))

    return mask


def rotate(x, y, x0, y0, angle):
    """ Rotate a line by ``angle`` [deg] around the point (x0, y0) """
    angle_rad = angle / 57.29578
    xnew = x0 + (x - x0) * np.cos(angle_rad) - (y - y0) * np.sin(angle_rad)
    ynew = y0 + (x - x0) * np.sin(angle_rad) + (y - y0) * np.cos(angle_rad)

    return xnew, ynew
