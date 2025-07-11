# -*- coding: utf-8 -*-
"""
Module contains the actual "Detector" effects.

As opposed to the Detector class and the DetectorManager (ex DetectorArray),
which are kept in the scopesim.detector subpackage. The effects here are used
to define the detector geometry and help creating the other classes mentioned
above. The effects here (or one of them) are what ScopeSim "actually needs" to
have defined in order to run.
"""

from typing import ClassVar

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

from ..optics.fov_volume_list import FovVolumeList
from .effects import Effect
from ..optics.image_plane_utils import (
    header_from_list_of_xy,
    calc_footprint,
    create_wcs_from_points,
)
from ..utils import (
    from_currsys,
    close_loop,
    figure_factory,
    quantity_from_table,
    unit_from_table,
    get_logger,
)

logger = get_logger(__name__)

__all__ = ["DetectorList", "DetectorWindow", "DetectorList3D"]


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
            filename : "detector_list.dat"
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

    dims = "xy"  # 2 spatial dimensions
    z_order: ClassVar[tuple[int, ...]] = (90, 290, 390, 490)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "pixel_scale": "!INST.pixel_scale",  # arcsec
            "active_detectors": "all",
        }
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
        if not isinstance(obj, FovVolumeList):
            return obj

        logger.debug("apply_to got kwargs: %s", kwargs)
        shrink_axis, shrink_values = self._get_fov_limits(
            pixel_scale=kwargs.get("pixel_scale", None))
        obj.shrink(axis=shrink_axis, values=shrink_values)

        return obj

    @property
    def pixel_size(self) -> u.Quantity[u.mm]:
        """Return size of one pixel in mm."""
        pixel_size = np.min(
            quantity_from_table("pixel_size", self.active_table, u.mm))
        return pixel_size

    @property
    def pixel_scale_mm(self) -> u.Equivalency:
        """Return pixel scale (mm / pix) as equivalency."""
        return u.pixel_scale(self.pixel_size / u.pixel)

    def _get_pixel_scale_arcsec(self, pixel_scale=None) -> u.Equivalency:
        """To allow overriding defaults, but use defaults in property."""
        pixel_scale = u.pixel_scale(
            from_currsys(
                pixel_scale if pixel_scale is not None else self.meta["pixel_scale"],
                self.cmds,
            ) << u.arcsec / u.pixel
        )
        return pixel_scale

    @property
    def pixel_scale_arcsec(self) -> u.Equivalency:
        """Return pixel scale (arcsec / pix) as equivalency."""
        return self._get_pixel_scale_arcsec()

    def _get_fov_limits(self, pixel_scale=None):
        mm_points = self._get_corner_points()[:, :2]
        pixel_scale = self._get_pixel_scale_arcsec(pixel_scale)

        # The line below combines the pixesl <=> arcsec scale with the
        # pixel <=> mm scale (both are astropy equivalencies, which can be
        # simply added together), to allow an indirect conversion from mm via
        # pixel to arcsec, converting the detector's mm footprint to an on-sky
        # FOV limit. Conceptually, this should involve the plate scale, but
        # ScopeSim is still somewhat inconsistent / redundant in defining
        # pixel size, pixel scale and plate scale.
        # TODO: Consider using plate scale here once that's generally been
        #       refactored to be more consistent everywhere.
        with u.set_enabled_equivalencies(pixel_scale + self.pixel_scale_mm):
            sky_points = (mm_points << u.pixel).to_value(u.arcsec)

        axis = ["x", "y"]
        values = (tuple(zip(sky_points.min(axis=0), sky_points.max(axis=0))))
        return axis, values

    def _get_corner_points(self):
        tbl = self.active_table
        with u.set_enabled_equivalencies(self.pixel_scale_mm):
            cens = {
                dim: (tbl[f"{dim}_cen"].data.astype(float) *
                      unit_from_table(f"{dim}_cen", tbl, u.mm) << u.mm)
                for dim in self.dims
            }

            ds = {
                dim: (0.5 * tbl[f"{dim}_size"].data.astype(float) *
                      unit_from_table(f"{dim}_size", tbl, u.mm) << u.mm)
                for dim in self.dims
            }

        det_mins = {dim: np.min(cens[dim] - ds[dim]) for dim in self.dims}
        det_maxs = {dim: np.max(cens[dim] + ds[dim]) for dim in self.dims}

        points = np.array([
            [det_mins[dim].to_value(u.mm) for dim in self.dims],
            [det_maxs[dim].to_value(u.mm) for dim in self.dims],
        ])

        return points * u.mm

    def fov_grid(self, which="edges", **kwargs):
        """Return an ApertureMask object. kwargs are "pixel_scale" [arcsec]."""
        # This has really been taken care of by apply_to now
        # TODO v1.0: finally rm this completely
        raise AttributeError("The DetectorList.fov_grid() method has been "
                             "removed in v0.10.0.")
        # aperture_mask = None
        # if which == "edges":
        #     self.meta.update(kwargs)
        #     self.meta = from_currsys(self.meta, self.cmds)

        #     hdr = self.image_plane_header
        #     xy_mm = calc_footprint(hdr, "D")
        #     pixel_size = hdr["CDELT1D"]              # mm
        #     pixel_scale = self.meta["pixel_scale"]   # ["]

        #     # x["] = x[mm] * ["] / [mm]
        #     xy_sky = xy_mm * pixel_scale / pixel_size

        #     aperture_mask = ApertureMask(array_dict={"x": xy_sky[:, 0],
        #                                              "y": xy_sky[:, 1]},
        #                                  pixel_scale=pixel_scale)

        # return aperture_mask

    @property
    def image_plane_header(self):
        """Create and return the Image Plane Header."""
        # FIXME: Heavy property.....
        # TODO: refactor this using (parts of) 3D implementation below
        tbl = self.active_table
        pixel_size = np.min(quantity_from_table("pixel_size", tbl, u.mm))
        x_size_unit = unit_from_table("x_size", tbl, u.mm)
        y_size_unit = unit_from_table("y_size", tbl, u.mm)
        # This is mm everywhere in the IRDB, but could be something else...
        x_cen_unit = unit_from_table("x_cen", tbl, u.mm)
        y_cen_unit = unit_from_table("y_cen", tbl, u.mm)

        xcen = tbl["x_cen"].data.astype(float) * x_cen_unit
        ycen = tbl["y_cen"].data.astype(float) * y_cen_unit
        dx = 0.5 * tbl["x_size"].data.astype(float) * x_size_unit
        dy = 0.5 * tbl["y_size"].data.astype(float) * y_size_unit

        pixel_scale = u.pixel_scale(pixel_size / u.pixel)
        with u.set_enabled_equivalencies(pixel_scale):
            xcen <<= u.mm
            ycen <<= u.mm
            dx <<= u.mm
            dy <<= u.mm

        x_det_min = np.min(xcen - dx)
        x_det_max = np.max(xcen + dx)
        y_det_min = np.min(ycen - dy)
        y_det_max = np.max(ycen + dy)

        x_det = [x_det_min.to_value(u.mm), x_det_max.to_value(u.mm)]
        y_det = [y_det_min.to_value(u.mm), y_det_max.to_value(u.mm)]

        pixel_size = pixel_size.to_value(u.mm)
        hdr = header_from_list_of_xy(x_det, y_det, pixel_size, "D")
        hdr["IMGPLANE"] = self.image_plane_id

        return hdr

    @property
    def active_table(self):
        """Create and return the active table."""
        # FIXME: Heavy property.....
        if self.meta["active_detectors"] == "all":
            tbl = self.table
        elif isinstance(self.meta["active_detectors"], (list, tuple)):
            mask = [det_id in self.meta["active_detectors"]
                    for det_id in self.table["id"]]
            tbl = self.table[mask]
        else:
            raise ValueError("Could not determine which detectors are active: "
                             f"{self.meta['active_detectors']}, {self.table},")
        tbl = from_currsys(tbl, self.cmds)

        return tbl

    def detector_headers(self, ids=None):
        """Create detector headers from active detectors or given IDs."""
        # Note: the ids argument is not used anywhere in ScopeSim.

        if ids is not None and all(isinstance(ii, int) for ii in ids):
            self.meta["active_detectors"] = list(ids)

        tbl = from_currsys(self.active_table, self.cmds)
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
        """Plot the detector layout."""
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

    @property
    def image_plane_id(self) -> int:
        """Get ID of the corresponding image plane."""
        return self.meta["image_plane_id"]


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

        array_dict = {
            "id": [0],
            "x_cen": [x],
            "y_cen": [y],
            "x_size": [width],
            "y_size": [height],
            "angle": [angle],
            "gain": [gain],
            "pixel_size": [pixel_size],
        }

        super().__init__(array_dict=array_dict, **params)


class DetectorList3D(DetectorList):
    """Pseudo-detector for simple IFU mode.

    .. versionadded:: 0.10.0
    """

    dims = "xyz"  # 2 spatial dimensions, 1 spectral dimension

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "pixel_scale": "!INST.pixel_scale",  # arcsec
            "active_detectors": "all",
            "wavelen": "!OBS.wavelen",
            "dwave": "!SIM.spectral.spectral_bin_width",
        }
        self.meta.update(params)
        self.meta.update(kwargs)

    @property
    def waverange(self) -> u.Quantity[u.um]:
        """Return wavelength range in um [wave_min, wave_max]."""
        # TODO: dynamic wave unit ?
        wave_points = self._get_corner_points()[:, 2]
        wave_mid = from_currsys(self.meta["wavelen"], self.cmds) * u.um
        dwave = from_currsys(self.meta["dwave"], self.cmds) * u.um / u.pixel

        with u.set_enabled_equivalencies(self.pixel_scale_mm):
            waverange = wave_points.to(u.pixel) * dwave + wave_mid
        return waverange

    def _get_fov_limits(self, pixel_scale=None):
        mm_points = self._get_corner_points()[:, :2]
        pixel_scale = self._get_pixel_scale_arcsec(pixel_scale)

        # See parent class for explanation of this operation.
        with u.set_enabled_equivalencies(pixel_scale + self.pixel_scale_mm):
            sky_points = (mm_points << u.pixel).to_value(u.arcsec)

        axis = ["x", "y", "wave"]
        values = (
            *zip(sky_points.min(axis=0), sky_points.max(axis=0)),
            tuple(self.waverange.to_value(u.um).round(7))
        )

        return axis, values  # [arcsec, arcsec, um]

    @property
    def image_plane_header(self):
        """Create and return the Image Plane Header."""
        # FIXME: Heavy property.....
        points = self._get_corner_points().value
        new_wcs, naxis = create_wcs_from_points(
            points, self.pixel_size.to_value(u.mm), "D")

        hdr = fits.Header()
        hdr["NAXIS"] = 3
        hdr["NAXIS1"] = int(naxis[0])
        hdr["NAXIS2"] = int(naxis[1])
        hdr["NAXIS3"] = int(naxis[2])
        hdr.update(new_wcs.to_header())
        hdr["IMGPLANE"] = self.image_plane_id

        return hdr

    def detector_headers(self, ids=None):
        """Override for simplified 3D case."""
        return [self.image_plane_header]
