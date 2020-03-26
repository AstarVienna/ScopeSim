import numpy as np
from astropy.io import fits
from astropy.table import Table

from .effects import Effect
from ..optics.spectral_trace import SpectralTrace
from ..optics import spectral_trace_utils as spt_utils
from ..utils import from_currsys, check_keys
from ..base_classes import FieldOfViewBase


class SpectralTraceListOld(Effect):
    """
    List of spectral trace geometries for the detector plane

    Should work in concert with an ApertureList (or ApertureMask) object and a
    DetectorList object

    Spectral trace patterns are to be kept in a ``fits.HDUList`` with one or
    more ``fits.BinTableHDU`` extensions, each one describing the geometry of a
    single trace. The 1st extension should be a ``BinTableHDU`` connecting the
    traces to the correct ``Aperture`` and ``ImagePlane`` objects.

    ..todo:: Add documentation describing the format of the fits objects

    The ``fits.HDUList`` objects can be loaded use one of the two keywords:
    - ``filename``: for on disk FITS files, or
    - ``hdulist``: for in-memory ``fits.HDUList`` objects

    """
    def __init__(self, **kwargs):
        super(SpectralTraceList, self).__init__(**kwargs)
        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        params = {"pixel_scale": "!INST.pixel_scale",       # [arcsec / pix]}
                  "plate_scale": "!INST.plate_scale",       # [arcsec / mm]
                  "wave_min":    "!SIM.spectral.wave_min",  # [um]
                  "wave_max":    "!SIM.spectral.wave_max",  # [um]
                  "x_colname":   "x",
                  "y_colname":   "y",
                  "s_colname":   "s",
                  "wave_colname": "wavelength",
                  "col_number_start": 0,
                  "dwave": 0.002,   # [um] for finding the best fit dispersion
                  "invalid_value": None}   # for dodgy trace file values

        self.meta["z_order"] = [70]
        self.meta.update(params)
        self.meta.update(kwargs)
        self.apply_to_classes = []

        # required_keys = ["pixel_scale", "plate_scale"]
        # check_keys(self.meta, required_keys, action="error")

    def fov_grid(self, which="waveset", **kwargs):
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta)

        if which == "waveset":
            return self.get_waveset()

        elif which == "edges":
            check_keys(kwargs, ["sky_header", "wave_min", "wave_max"], "error")
            return self.get_fov_headers(kwargs["sky_header"])

    def get_waveset(self):
        self.meta = from_currsys(self.meta)
        params = {}
        params.update(self.meta)

        ext_data = self._file[0].header["EDATA"]
        tbls = [hdu.data for hdu in self._file[ext_data:]]  # all trace tables

        if "wave_min" not in params or "wave_max" not in params:
            wave_colname = self.meta["wave_colname"]
            wave_min = np.min([tbl[wave_colname].max() for tbl in tbls])
            wave_max = np.max([tbl[wave_colname].max() for tbl in tbls])
            params.update({"wave_min": wave_min, "wave_max": wave_max})
        else:
            wave_min, wave_max = params["wave_min"], params["wave_max"]

        disp, waverange = spt_utils.get_max_dispersion(tbls, **params)
        um_per_pix = self.meta["pixel_scale"] / self.meta["plate_scale"] / disp
        waveset = spt_utils.pixel_wavelength_edges(um_per_pix, waverange,
                                                   wave_min, wave_max)
        return waveset

    def get_fov_headers(self, sky_header):
        self.meta = from_currsys(self.meta)

        dwave = self.meta["dwave"]
        pixel_size = self.meta["pixel_scale"] / self.meta["plate_scale"]

        # based on on-sky aperture_id, get relevant traces and image_plane_ids
        cat_ext = self._file[0].header["ECAT"]
        cat_tbl = self._file[cat_ext].data
        mask = cat_tbl["aperture_id"] == sky_header["ID"]
        ext_ids = cat_tbl["extension_id"][mask]
        implane_ids = cat_tbl["image_plane_id"][mask]

        wave_min = self.meta["wave_min"]
        wave_max = self.meta["wave_max"]
        # for each trace get the max_dispersion waveset
        headers_list = []
        for ext, implane in zip(ext_ids, implane_ids):
            trace = SpectralTrace(self._file[ext].data, **self.meta)

            if trace.wave_min < wave_max and trace.wave_max > wave_min:
                trace.get_trace_curves(pixel_size=pixel_size,
                                       wave_min=wave_min,
                                       wave_max=wave_max)
                hdrs = trace.get_curve_headers(pixel_size=pixel_size)

                # combine sky header with monochromatic trace curve header
                for mtc_hdr in hdrs:
                    mtc_hdr["IMGPLANE"] = implane
                    mtc_plate_scale = mtc_hdr["PLATESCL"] / 3600.  # [arcsec]
                    pixel_D_scale = sky_header["CDELT1"] / mtc_plate_scale
                    mtc_hdr["CDELT1D"] = pixel_D_scale
                    mtc_hdr["CDELT2D"] = pixel_D_scale
                    mtc_hdr["CRPIX2D"] = sky_header["NAXIS2"] * 0.5
                    # ..todo:: assumption here is that they are on the same pixel scale - bad assumption!
                    mtc_hdr["NAXIS2"] = sky_header["NAXIS2"]
                    mtc_hdr.update(sky_header)
                    # ??? speed up potential ???
                    # mtc_hdr = combine_sky_slit_headers(sky_header, mtc_hdr)

                headers_list += hdrs

        return headers_list

    def apply_to(self, fov):
        if isinstance(fov, self.apply_to_classes):
            pass

        return fov

    def plot(self, wave_min=None, wave_max=None, **kwargs):
        import matplotlib.pyplot as plt

        if wave_min is None:
            wave_min = self.meta["wave_min"]
        if wave_max is None:
            wave_max = self.meta["wave_max"]

        ext_data = self._file[0].header["EDATA"]
        for ext, c in zip(self._file[ext_data:], "rgbymc" * len(self._file)):
            tbl = Table(ext.data)
            n_traces = len([col for col in tbl.colnames
                            if self.meta["y_colname"] in col])
            col_start = self.meta["col_number_start"]
            for i in range(col_start, col_start + n_traces):
                wave = tbl[self.meta["wave_colname"]]
                x = tbl[self.meta["x_colname"] + str(i)]
                y = tbl[self.meta["y_colname"] + str(i)]
                if np.min(wave) > wave_max or np.max(wave) < wave_min:
                    continue

                mask = (wave >= wave_min) * (wave <= wave_max)
                plt.plot(x[mask], y[mask], c, alpha=0.5)
