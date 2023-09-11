import numpy as np

from ..base_classes import PoorMansHeader


class MonochromeTraceCurve:
    """
    Contains coordinates along a monochromatic section of a spectral trace.

    Used to generate detector plane position WCS info for FOVs in spectral mode
    """

    def __init__(self, x, y, s, wave_min, wave_max, **kwargs):
        self.wave_min = wave_min
        self.wave_max = wave_max
        self.x = x
        self.y = y
        self.s = s

        self.meta = {"rotation": 0,
                     "shear": 0}
        self.meta.update(kwargs)

        self._header = None

    @property
    def header(self):
        if "pixel_size" in self.meta:
            hdr = self.get_header(self.meta["pixel_size"])
        else:
            hdr = self._header

        return hdr

    def get_header(self, pixel_size):

        dx = self.x[-1] - self.x[0]     # [mm]
        dy = self.y[-1] - self.y[0]     # [mm]
        ds = self.s[-1] - self.s[0]     # [mm]
        len_mm = (dx**2 + dy**2)**0.5   # [mm]
        arcsec_per_mm = ds / len_mm     # [arcsec / mm]

        x_sign = 1 if dx * ds >= 0 else -1
        y_sign = 1 if dy * ds >= 0 else -1

        hdr_dict = {"NAXIS": 2,
                    "NAXIS1": int(round(len_mm / pixel_size)),
                    "NAXIS2": 0,
                    "CTYPE1D": "LINEAR",
                    "CTYPE2D": "LINEAR",
                    "CUNIT1D": "mm",
                    "CUNIT2D": "mm",
                    "CDELT1D": pixel_size * x_sign,    # [mm / pix]
                    "CDELT2D": pixel_size * y_sign,    # [mm / pix]
                    "CRVAL1D": self.x[0],
                    "CRVAL2D": self.y[0],
                    "CRPIX1D": 0,
                    "CRPIX2D": 0,
                    }

        rad2deg = np.pi / 180.
        c = np.cos(self.meta["rotation"] * rad2deg)
        s = np.sin(self.meta["rotation"] * rad2deg)
        t = np.tan(self.meta["shear"] * rad2deg)
        pc_dict = {"PC1_1D": c, "PC1_2D": s + c*t,
                   "PC2_1D": -s, "PC2_2D": c + s*t}
        hdr_dict.update(pc_dict)

        # ..todo: deal with these aspects
        # curvature

        # self._header = fits.Header(hdr_dict)
        self._header = PoorMansHeader(hdr_dict)
        self._header["WAVE_MIN"] = self.wave_min
        self._header["WAVE_MID"] = 0.5 * (self.wave_min + self.wave_max)
        self._header["WAVE_MAX"] = self.wave_max
        self._header["PLATESCL"] = (arcsec_per_mm,
                                    "[arcsec mm-1] Along slit")
        self._header["ROTANGD"] = (self.meta["rotation"] * rad2deg,
                                   "[deg] w.r.t. focal plane")
        self._header["SKEWANGD"] = (self.meta["shear"] * rad2deg,
                                    "[deg] w.r.t. trace trajectory")

        return self._header
