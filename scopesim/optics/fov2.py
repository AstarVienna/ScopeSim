import numpy as np

from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

from . import fov_utils as fu
from . import image_plane_utils as imp_utils

from ..base_classes import SourceBase, FieldOfViewBase, PoorMansHeader
from .. import utils


class FieldOfView(FieldOfViewBase):
    """
    A FOV is spectro-spatial volume cut out of a Source object.

    Flux units after extracting the fields from the Source are in ph/s/pixel

    The initial header should contain an on-sky WCS description:
    - CDELT-, CUNIT-, NAXIS- : for pixel scale and size (assumed CUNIT in deg)
    - CRVAL-, CRPIX- : for positioning the final image
    - CTYPE- : is assumed to be "RA---TAN", "DEC---TAN"

    and an image-plane WCS description
    - CDELT-D, CUNIT-D, NAXISn : for pixel scale and size (assumed CUNIT in mm)
    - CRVAL-D, CRPIX-D : for positioning the final image
    - CTYPE-D : is assumed to be "LINEAR", "LINEAR"

    The wavelength range is given by waverange

    """

    def __init__(self, header, waverange, **kwargs):
        self.meta = {"id": None,
                     "wave_min": utils.quantify(waverange[0], u.um),
                     "wave_max": utils.quantify(waverange[1], u.um),
                     "area": 0 * u.m**2,
                     "sub_pixel": "!SIM.sub_pixel.flag",
                     "distortion": {"scale": [1, 1],
                                    "offset": [0, 0],
                                    "shear": [1, 1],
                                    "rotation": 0,
                                    "radius_of_curvature": None},
                     "conserve_image": True,
                     "trace_id": None,
                     "aperture_id": None,
                     }
        self.meta.update(kwargs)

        if not any([utils.has_needed_keywords(header, s) for s in ["", "S"]]):
            raise ValueError("header must contain a valid sky-plane WCS: {}"
                             "".format(dict(header)))
        if not utils.has_needed_keywords(header, "D"):
            raise ValueError("header must contain a valid image-plane WCS: {}"
                             "".format(dict(header)))

        if isinstance(header, PoorMansHeader):
            self.hdu = fits.ImageHDU()
            self.hdu.header.update(header)
        else:
            self.hdu = fits.ImageHDU(header=header)
        self.hdu.header["NAXIS"] = 2
        self.hdu.header["NAXIS1"] = header["NAXIS1"]
        self.hdu.header["NAXIS2"] = header["NAXIS2"]

        self.image_plane_id = 0
        self.fields = []
        self.spectra = []

        self.cube = None        # IFU, long-lit, Slicer-MOS
        self.image = None       # Imagers
        self.spectrum = None    # Fibre-fed MOS

        self._waverange = None
        self._wavelength = None
        self._volume = None

    def extract_from(self, src):
        """ ..assumption: Bandpass has been applied"""

        if not isinstance(src, SourceBase):
            raise ValueError("source must be a Source object: {}"
                             "".format(type(src)))

        fields_in_fov = [field for field in src.fields
                         if fu.is_field_in_fov(self.hdu.header, field)]

        spec_refs = []
        volume = self.volume()
        for i in range(len(fields_in_fov)):
            fld = fields_in_fov[i]
            if isinstance(fld, Table):
                fields_in_fov[i] = fu.extract_area_from_table(fld, volume)
                spec_refs += list(np.unique(fields_in_fov[i] ["ref"]))

            elif isinstance(fld, fits.ImageHDU) and fld.header["NAXIS"] == 2:
                fields_in_fov[i] = fu.extract_area_from_imagehdu(fld, volume)
                ref = fld.header.get("SPEC_REF")
                if ref is not None:
                    spec_refs += [ref]

            elif isinstance(fld, fits.ImageHDU) and fld.header["NAXIS"] == 3:
                fields_in_fov[i] = fu.extract_area_from_imagehdu(fld, volume)

        waves = volume["waves"] * u.Unit(volume["wave_unit"])
        spectra = {ref: fu.extract_range_from_spectrum(src.spectra[ref], waves)
                   for ref in np.unique(spec_refs)}

        self.fields = fields_in_fov
        self.spectra = spectra

    def view(self, sub_pixel=None):
        return None

    def make_spectrum(self):
        # This is needed for when we do incoherent MOS instruments.
        # Each fibre doesn't care about the spatial information.
        return None

    def make_image(self):
        # Used for imaging
        return None

    def make_cube(self):
        # Used for IFUs, slit spectrographs, and coherent MOSs (e.g.KMOS)
        return None

    def volume(self, wcs_prefix=""):
        xs, ys = imp_utils.calc_footprint(self.header, wcs_suffix=wcs_prefix)
        wave_corners = self.waverange
        self._volume = {"xs": [min(xs), max(xs)],
                        "ys": [min(ys), max(ys)],
                        "waves": self.waverange,
                        "xy_unit": "mm" if wcs_prefix == "D" else "deg",
                        "wave_unit": "um"}
        return self._volume

    @property
    def header(self):
        return self.hdu.header

    @property
    def data(self):
        if self.hdu.data is None:
            self.view(self.meta["sub_pixel"])
        return self.hdu.data

    @property
    def corners(self):
        sky_corners = imp_utils.calc_footprint(self.header)
        imp_corners = imp_utils.calc_footprint(self.header, "D")
        return sky_corners, imp_corners

    @property
    def waverange(self):
        """Returns wavelength range in um [wave_min, wave_max]"""
        if self._waverange is None:
            wave_min = utils.quantify(self.meta["wave_min"], u.um).value
            wave_max = utils.quantify(self.meta["wave_max"], u.um).value
            self._waverange = [wave_min, wave_max]
        return self._waverange

    @property
    def wavelength(self):
        """Returns central wavelength in um"""
        if self._wavelength is None:
            self._wavelength = np.average(self.waverange)
        return self._wavelength

    def __repr__(self):
        msg = "FOV id: {}, with dimensions ({}, {})\n" \
              "".format(self.meta["id"], self.header["NAXIS1"],
                        self.header["NAXIS2"])
        msg += "Sky centre: ({}, {})\n" \
               "".format(self.header["CRVAL1"], self.header["CRVAL2"])
        msg += "Image centre: ({}, {})\n" \
               "".format(self.header["CRVAL1D"], self.header["CRVAL2D"])
        msg += "Wavelength range: ({}, {})um\n" \
               "".format(self.meta["wave_min"], self.meta["wave_max"])

        return msg
