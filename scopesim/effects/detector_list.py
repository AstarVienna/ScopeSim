"""TBA."""

import logging

import numpy as np
from astropy import units as u
from astropy.table import Table

from ..base_classes import FOVSetupBase
from .effects import Effect
from .apertures import ApertureMask
from .. import utils
from ..utils import close_loop, figure_factory
from ..optics.image_plane_utils import header_from_list_of_xy, calc_footprint

__all__ = ["DetectorList", "DetectorWindow"]


class DetectorList(Effect):
    """
    A description of detector positions and properties.

    The list of detectors must have the following table columns
    ::

        id   x_cen   y_cen  x_size  y_size  pixel_size  angle    gain

    where:

    * "id" is a reference id for the chip (fits header EXTNAME)
    * "x_cen" and "y_cen" [mm] are the physical coordinates of centre of the chip on the detector plane
    * "x_size", "y_size" [mm, pixel] are the width/height of the chip
    * "pixel_size" [mm] is the physical size of pixels in the detector
    * "angle" [deg] is the rotation of the detector relative to the x-axis
    * "gain" [e-/ADU] is the conversion factor for electrons (photons) to ADUs

    The units for each column (except ``id``) must be given in the meta data
    using the format ``<colname>_unit``. E.g. ``x_size_unit``.
    See examples below.

    .. note::
       Currently only the units specified below are accepted.

       For ``x(y)_size_unit``, acceptable units are ``mm``, ``pixel``

    Parameters
    ----------
    filename : str, optional
        Filename of the ASCII file with the detector description. See examples

    array_dict : dict
        Dict containing the detector description. See examples

    image_plane_id : int
        Which image plane the detector will look at (generally 0)

    Examples
    --------
    With the ``array_dict`` feature
    ::

        -   name: single_detector
            class: DetectorList
            kwargs:
                image_plane_id : 0
                array_dict:
                    id: [1]
                    x_cen: [0.]
                    y_cen: [0.]
                    x_size: [5.12]
                    y_size: [5.12]
                    pixel_size: [0.01]
                    angle: [0.]
                    gain: [1.0]
                x_cen_unit: mm
                y_cen_unit: mm
                x_size_unit: mm
                y_size_unit: mm
                pixel_size_unit: mm
                angle_unit: deg
                gain_unit: electron/adu


    Or referring to a table contained in a seperate ASCII file
    ::

        - name : full_detector_array
          class : DetectorList
          kwargs :
            filename : "detecotr_list.dat"
            active_detectors : [1, 3]
            image_plane_id : 0

    where the file detector_list.dat contains the following information
    ::

        # x_cen_unit : mm
        # y_cen_unit : mm
        # x_size_unit : pix
        # y_size_unit : pix
        # pixel_size_unit : mm
        # angle_unit : deg
        # gain_unit : electron/adu
        #
        id   x_cen   y_cen  x_size  y_size  pixel_size  angle    gain
        1   -63.94    0.00    4096    4096       0.015    0.0     1.0
        2     0.00    0.00    4096    4096       0.015   90.0     1.0
        3    63.94    0.00    4096    4096       0.015  180.0     1.0

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {"z_order": [90, 290, 390, 490],
                  "pixel_scale": "!INST.pixel_scale",      # arcsec
                  "active_detectors": "all",
                  "report_plot_include": True,
                  "report_table_include": True}
        self.meta.update(params)
        self.meta.update(kwargs)

        # for backwards compatibility
        new_colnames = {"xhw": "x_size",
                        "yhw": "y_size",
                        "pixsize": "pixel_size"}
        mult_cols = {"xhw": 2., "yhw": 2., "pixsize": 1.}
        if isinstance(self.table, Table):
            for col, new_name in new_colnames.items():
                if col in self.table.colnames:
                    self.table[col] = self.table[col] * mult_cols[col]
                    self.table.rename_column(col, new_name)
        if "x_size_unit" not in self.meta and "xhw_unit" in self.meta:
            self.meta["x_size_unit"] = self.meta["xhw_unit"]
        if "y_size_unit" not in self.meta and "yhw_unit" in self.meta:
            self.meta["y_size_unit"] = self.meta["yhw_unit"]

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, FOVSetupBase):

            hdr = self.image_plane_header
            xy_mm = calc_footprint(hdr, "D")
            pixel_size = hdr["CDELT1D"]              # mm
            pixel_scale = kwargs.get("pixel_scale", self.meta["pixel_scale"])   # ["]
            pixel_scale = utils.from_currsys(pixel_scale)

            # x["] = x[mm] * ["] / [mm]
            xy_sky = xy_mm * pixel_scale / pixel_size

            obj.shrink(axis=["x", "y"],
                       values=(tuple(zip(xy_sky.min(axis=0),
                                         xy_sky.max(axis=0)))))

            lims = np.array((xy_mm.min(axis=0), xy_mm.max(axis=0)))
            keys = ["xd_min", "xd_max", "yd_min", "yd_max"]
            obj.detector_limits = dict(zip(keys, lims.T.flatten()))

        return obj

    def fov_grid(self, which="edges", **kwargs):
        """Return an ApertureMask object. kwargs are "pixel_scale" [arcsec]."""
        logging.warning("DetectorList.fov_grid will be depreciated in v1.0")
        aperture_mask = None
        if which == "edges":
            self.meta.update(kwargs)
            self.meta = utils.from_currsys(self.meta)

            hdr = self.image_plane_header
            xy_mm = calc_footprint(hdr, "D")
            pixel_size = hdr["CDELT1D"]              # mm
            pixel_scale = self.meta["pixel_scale"]   # ["]

            # x["] = x[mm] * ["] / [mm]
            xy_sky = xy_mm * pixel_scale / pixel_size

            aperture_mask = ApertureMask(array_dict={"x": xy_sky[:, 0],
                                                     "y": xy_sky[:, 1]},
                                         pixel_scale=pixel_scale)

        return aperture_mask

    @property
    def image_plane_header(self):
        tbl = self.active_table
        pixel_size = np.min(utils.quantity_from_table("pixel_size", tbl, u.mm))
        x_unit = utils.unit_from_table("x_size", tbl, u.mm)
        y_unit = utils.unit_from_table("y_size", tbl, u.mm)

        xcen = tbl["x_cen"].data.astype(float)
        ycen = tbl["y_cen"].data.astype(float)
        dx = 0.5 * tbl["x_size"].data.astype(float)
        dy = 0.5 * tbl["y_size"].data.astype(float)

        scale_unit = 1        # either unitless to retain
        if "pix" in x_unit.name:
            scale_unit = u.mm / u.pix
            dx *= pixel_size.value
            dy *= pixel_size.value

        x_det_min = np.min(xcen - dx) * x_unit * scale_unit
        x_det_max = np.max(xcen + dx) * x_unit * scale_unit
        y_det_min = np.min(ycen - dy) * y_unit * scale_unit
        y_det_max = np.max(ycen + dy) * y_unit * scale_unit

        x_det = [x_det_min.to(u.mm).value, x_det_max.to(u.mm).value]
        y_det = [y_det_min.to(u.mm).value, y_det_max.to(u.mm).value]

        pixel_size = pixel_size.to(u.mm).value
        hdr = header_from_list_of_xy(x_det, y_det, pixel_size, "D")
        hdr["IMGPLANE"] = self.meta["image_plane_id"]

        return hdr

    @property
    def active_table(self):
        if self.meta["active_detectors"] == "all":
            tbl = self.table
        elif isinstance(self.meta["active_detectors"], (list, tuple)):
            mask = [det_id in self.meta["active_detectors"]
                    for det_id in self.table["id"]]
            tbl = self.table[mask]
        else:
            raise ValueError("Could not determine which detectors are active: "
                             f"{self.meta['active_detectors']}, {self.table},")
        tbl = utils.from_currsys(tbl)

        return tbl

    def detector_headers(self, ids=None):
        if ids is not None and all(isinstance(ii, int) for ii in ids):
            self.meta["active_detectors"] = list(ids)

        tbl = utils.from_currsys(self.active_table)
        hdrs = []
        for row in tbl:
            pixel_size = row["pixel_size"]
            xcen, ycen = float(row["x_cen"]), float(row["y_cen"])
            dx, dy = 0.5 * float(row["x_size"]), 0.5 * float(row["y_size"])

            if "pix" in self.meta["x_cen_unit"]:
                xcen, ycen = xcen * pixel_size, ycen * pixel_size
            if "pix" in self.meta["x_size_unit"]:
                dx, dy = dx * pixel_size, dy * pixel_size

            hdr = header_from_list_of_xy([xcen-dx, xcen+dx],
                                         [ycen-dy, ycen+dy],
                                         pixel_scale=pixel_size,
                                         wcs_suffix="D")
            if abs(row["angle"]) > 1E-4:
                sang = np.sin(row["angle"] / 57.29578)
                cang = np.cos(row["angle"] / 57.29578)
                hdr["PC1_1"] = cang
                hdr["PC1_2"] = sang
                hdr["PC2_1"] = -sang
                hdr["PC2_2"] = cang

            # hdr["GAIN"] = row["gain"]
            if "id" in row:
                hdr["DET_ID"] = row["id"]
                hdr["EXTNAME"] = f"DET_{row['id']}"

            row_dict = {col: row[col] for col in row.colnames}
            hdr.update(row_dict)
            hdrs.append(hdr)

        return hdrs

    def plot(self, axes=None):
        if axes is None:
            _, axes = figure_factory()

        for hdr in self.detector_headers():
            xy_mm = calc_footprint(hdr, "D")
            outline = np.array(list(close_loop(xy_mm)))
            axes.plot(outline[:, 0], outline[:, 1])
            axes.text(*xy_mm.mean(axis=0), hdr["ID"],
                      ha="center", va="center")

        axes.set_aspect("equal")
        axes.set_xlabel("Size [mm]")
        axes.set_ylabel("Size [mm]")

        return axes


class DetectorWindow(DetectorList):
    """
    For when a full DetectorList if too cumbersome.

    Parameters
    ----------
    pixel_size : float
        [mm pixel-1] Physical pixel size
    x, y : float
        [mm] Position of window centre relative to optical axis
    width, height=None : float
        [mm] Dimensions of window. If height is None, height=width
    angle : float, optional
        [deg] Rotation of window
    gain : float, optional
        [ADU/e-]
    units : str, optional
        [mm, pixel] Default "mm". Sets the input parameter units.
        If ``"pixel"``, (``x``, ``y``, ``width``, ``height``) are multiplied
        by ``pixel_size``

    """

    def __init__(self, pixel_size, x, y, width, height=None, angle=0, gain=1,
                 units="mm", **kwargs):

        if height is None:
            height = width

        params = {"orig_units": units,
                  "x_cen_unit": units,
                  "y_cen_unit": units,
                  "x_size_unit": units,
                  "y_size_unit": units,
                  "pixel_size_unit": "mm",
                  "angle_unit": "deg",
                  "gain_unit": "electron/adu",
                  "image_plane_id": 0}
        params.update(kwargs)

        tbl = Table(data=[[0], [x], [y], [width], [height],
                          [angle], [gain], [pixel_size]],
                    names=["id", "x_cen", "y_cen", "x_size", "y_size",
                           "angle", "gain", "pixel_size"])
        tbl.meta.update(params)

        super().__init__(table=tbl, **params)
