import numpy as np

from astropy.io import fits
from astropy.table import Table

from .effects import Effect
from ..optics.spectral_trace import SpectralTrace
from ..utils import from_currsys, check_keys
from ..optics.image_plane_utils import header_from_list_of_xy


class SpectralTraceList(Effect):
    """
    List of spectral trace geometries for the detector plane

    Should work in concert with an ApertureList (or ApertureMask) object and a
    DetectorList object

    Spectral trace patterns are to be kept in a ``fits.HDUList`` with one or
    more ``fits.BinTableHDU`` extensions, each one describing the geometry of a
    single trace. The 1st extension should be a ``BinTableHDU`` connecting the
    traces to the correct ``Aperture`` and ``ImagePlane`` objects.

    ..todo:: Add documentation describing the format of the fits objects

    The ``fits.HDUList`` objects can be loaded using one of these two keywords:

    - ``filename``: for on disk FITS files, or
    - ``hdulist``: for in-memory ``fits.HDUList`` objects

    The format and contents of the extensions in the HDUList (FITS file) object
    is listed below

    Input Data Format
    -----------------
    A trace list FITS file needs the following extensions:

    - 0 : PrimaryHDU [header]
    - 1 : BinTableHDU [header, data] : Overview table of all traces
    - 2..N : BinTableHDU [header, data] : Trace tables. One per spectral trace

    EXT 0 : PrimaryHDU
    ++++++++++++++++++
    Required Header Keywords:

    - ECAT : int : Extension number of overview table. Normally 1
    - EDATA : int : Extension number of first Trace table. Normally 2

    No data is required in this extension

    EXT 1 : BinTableHDU : Overview of traces
    ++++++++++++++++++++++++++++++++++++++++
    No special header keywords are required in this extension

    Required Table columns:

    - description : str : description of each each trace
    - extension_id : int : which extension is each trace in
    - aperture_id : int : which aperture matches this trace (e.g. MOS / IFU)
    - image_plane_id : int : on which image plane is this trace projected

    EXT 2 : BinTableHDU : Individual traces
    +++++++++++++++++++++++++++++++++++++++
    No special header keywords are required in this extension

    Required Table columns:

    - wavelength : float : [um] : wavelength of monochromatic aperture image
    - s0 .. s0 : float : [mm] : position along aperture perpendicular to trace
    - x0 .. xN : float : [mm] : x position of aperture image on focal plane
    - y0 .. yN : float : [arcsec] : y position of aperture image on focal plane



    """
    def __init__(self, **kwargs):
        super(SpectralTraceList, self).__init__(**kwargs)

        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        self.meta.update({"pixel_scale": "!INST.pixel_scale",  # [arcsec / pix]}
                          "plate_scale": "!INST.plate_scale",  # [arcsec / mm]
                          "wave_min": "!SIM.spectral.wave_min",  # [um]
                          "wave_max": "!SIM.spectral.wave_max",  # [um]
                          "x_colname": "x",
                          "y_colname": "y",
                          "s_colname": "s",
                          "wave_colname": "wavelength",
                          "col_number_start": 0,
                          "dwave": 0.002,  # [um] for finding the best fit dispersion
                          "invalid_value": None  # for dodgy trace file values
                          })
        self.meta.update(kwargs)
        self.meta["z_order"] = [70, 270]

        if self._file is not None:
            self.ext_data = self._file[0].header["EDATA"]
            self.ext_cat = self._file[0].header["ECAT"]
            self.catalog = Table(self._file[self.ext_cat].data)
            self.spectral_traces = self.make_spectral_traces()

    def make_spectral_traces(self):
        spec_traces = []
        for row in self.catalog:
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            hdu = self._file[row["extension_id"]]
            spec_traces += [SpectralTrace(hdu, **params)]

        return spec_traces

    def fov_grid(self, which="waveset", **kwargs):
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta)

        if which == "waveset":
            if "pixel_size" not in kwargs:
                kwargs["pixel_size"] = None
            return self.get_waveset(kwargs["pixel_size"])

        elif which == "edges":
            check_keys(kwargs, ["sky_header", "det_header",
                                "wave_min", "wave_max",
                                "pixel_scale", "plate_scale"], "error")
            return self.get_fov_headers(**kwargs)

    def get_waveset(self, pixel_size=None):
        if pixel_size is None:
            pixel_size = self.meta["pixel_scale"] / self.meta["plate_scale"]

        wavesets = [spt.get_pixel_wavelength_edges(pixel_size)
                    for spt in self.spectral_traces]

        return wavesets

    def get_fov_headers(self, sky_header, **kwargs):
        fov_headers = []
        for spt in self.spectral_traces:
            fov_headers += spt.fov_headers(sky_header=sky_header, **kwargs)

        return fov_headers

    def apply_to(self, fov):
        return fov

    @property
    def footprint(self):
        xs, ys = [], []
        for spt in self.spectral_traces:
            xi, yi = spt.footprint
            xs += xi
            ys += yi

        xs = [np.min(xs), np.max(xs), np.max(xs), np.min(xs)]
        ys = [np.min(ys), np.min(ys), np.max(ys), np.max(ys)]

        return xs, ys

    @property
    def image_plane_header(self):
        x, y = self.footprint
        pixel_scale = from_currsys(self.meta["pixel_scale"])
        hdr = header_from_list_of_xy(x, y, pixel_scale, "D")

        return hdr

    def plot(self, wave_min=None, wave_max=None, **kwargs):
        if self.spectral_traces is not None:
            clrs = "rgbcymk" * (1 + len(self.spectral_traces) // 7)
            for spt, c in zip(self.spectral_traces, clrs):
                spt.plot(wave_min, wave_max, c=c)

    def __repr__(self):
        return "\n".join([spt.__repr__() for spt in self.spectral_traces])
