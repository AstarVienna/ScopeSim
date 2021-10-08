import numpy as np
from scipy.interpolate import interp1d

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
                     "wave_bin_n": 1,
                     "wave_bin_type": "linear",
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
        self.spectra = {}

        self.cube = None        # IFU, long-lit, Slicer-MOS
        self.image = None       # Imagers
        self.spectrum = None    # Fibre-fed MOS

        self._waverange = None
        self._wavelength = None
        self._volume = None

    def extract_from(self, src):
        """ ..assumption: Bandpass has been applied

        .. note:: Spectra are cut and copied from the original Source object
            They are in original units. ph/s/pix comes in the make_**** methods

        """

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
        """
        Used for imaging

        Image units are ph s-1 m-2

        .. note:: ``self.make_image()`` does NOT store anything in ``self.image``

            See make_cube for an explanation

        1. Make waveset and canvas image
            make canvas image from NAXIS1,2 from fov.header

        2. Find Cube fields
            collapse cube along wavelength axis
            rescale and reorient image
            add cube image to canvas image

        3. Find Image fields
            rescale and reorient images
            sum spectra between wavelength edges
            multiply image by summed spectrum
            add image to canvas image

        4. Find Table fields
            sum spectra between wavelength edges
            add summed flux at x,y position in canvas image

        """

        if self.cube is not None:
            image = cube.sum(axis=0)
            canvas_image_hdu = fits.ImageHDU(data=image,
                                             header=self.hdu.header)
            return canvas_image_hdu

        # 1. Make waveset and canvas image
        fov_waveset = self.waveset
        naxis1, naxis2 = self.hdu.header["NAXIS1"], self.hdu.header["NAXIS2"]
        canvas_image_hdu = fits.ImageHDU(data=np.zeros((naxis2, naxis1)),
                                         header=self.hdu.header)

        # 2. Find Cube fields
        for field in self.cube_fields:
            field.data.sum(axis=0)

            1 / 0

        return None

    def make_cube(self):
        """
        Used for IFUs, slit spectrographs, and coherent MOSs (e.g.KMOS)

        Returned cube units are ph s-1 voxel-1

        .. note:: ``self.make_cube()`` does NOT store anything in ``self.cube``

            self.cube and self.make_cube() are deliberately keep seperately
            so that self.cube will not be accidently overwritten by a rogue
            call from an Effect object.

            All Effect objects should specifically test whether
            ``self.cube is None`` before assigning a new cube it

        The cube is made with these steps:

        1. Make waveset and canvas cube
            if at least one cube:
                set waveset to equal largest cube waveset
            else:
                make waveset from self.meta values
            make canvas cube based on waveset of largest cube and NAXIS1,2 from fov.header

        2. Find Cube fields
            rescale and reorient cubes
            interp1d smaller cubes with waveset
            add cubes to cavas cube

        3. Find Image fields
            rescale and reorient images
            evaluate spectra at waveset
            expand image by spectra to 3D form
            add image cubes to canvas cube

        4. Find Table fields
            evaluate spectra at waveset
            add spectrum at x,y position in canvas cube

        PHOTLAM = ph/s/m2/um
        original source fields are in units of:
        - tables: PHOTLAM (spectra)
        - images: PHOTLAM (spectra)
        - cubes: PHOTLAM * bin_width

        .. warning:: Images and Cubes need to be in units of pixel-1 (or voxel-1), not arcsec-2

        """
        # 1. Make waveset and canvas cube
        fov_waveset = self.waveset
        specs = {ref: spec(fov_waveset)     # * bin_widths * area
                 for ref, spec in self.spectra.items()}

        # make canvas cube based on waveset of largest cube and NAXIS1,2 from fov.header
        naxis1, naxis2 = self.hdu.header["NAXIS1"], self.hdu.header["NAXIS2"]
        naxis3 = len(fov_waveset)
        canvas_cube_hdu = fits.ImageHDU(data=np.zeros((naxis3, naxis2, naxis1)),
                                        header=self.hdu.header)
        canvas_cube_hdu.header["BUNIT"] = "ph s-1 cm-2 AA-1"

        # 2. Add Cube fields
        for field in self.cube_fields:
            field_waveset = fu.get_cube_waveset(field.header,
                                                return_quantity=True)
            field_interp = interp1d(field_waveset.to(u.um).value,
                                    field.data, axis=0, kind="linear")
            field_data = field_interp(fov_waveset.value)
            field_unit = field.header.get("BUNIT", "ph s-1 cm-2 AA-1")
            flux_scale_factor = u.Unit(field_unit).to("ph s-1 cm-2 AA-1")
            field_hdu = fits.ImageHDU(data=field_data * flux_scale_factor,
                                      header=field.header)
            canvas_cube_hdu = imp_utils.add_imagehdu_to_imagehdu(field_hdu,
                                                                 canvas_cube_hdu)

        # 3. Find Image fields
        for field in self.image_fields:
            canvas_image_hdu = fits.ImageHDU(data=np.zeros((naxis2, naxis1)),
                                             header=self.hdu.header)
            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(field,
                                                               canvas_image_hdu)
            spec = specs[field.header["SPEC_REF"]]
            field_cube = canvas_image_hdu.data[None, :, :] * spec[:, None, None]  # 2D * 1D -> 3D
            canvas_cube_hdu.data += field_cube.value

        # 4. Find Table fields
        for field in self.table_fields:
            for xsky, ysky, ref, weight in field:  # field is a Table and the iterable is the row vector
                # x, y are ALWAYS in arcsec - crval is in deg
                xpix, ypix = imp_utils.val2pix(self.hdu.header, xsky / 3600, ysky / 3600)
                if utils.from_currsys(self.meta["sub_pixel"]):
                    xs, ys, fracs = imp_utils.sub_pixel_fractions(xpix, ypix)
                    for i, j, k in zip(xs, ys, fracs):
                        wk = weight * k
                        canvas_cube_hdu.data[:, j, i] += specs[ref].value * wk
                else:
                    x, y = int(xpix), int(ypix)
                    flux_vector = specs[ref].value * weight
                    canvas_cube_hdu.data[:, y, x] += flux_vector

        # 5. Convert from PHOTLAM to ph/s/voxel
        #    PHOTLAM = ph/s/cm-2/AA
        #    area = m2, fov_waveset = um

        area = utils.from_currsys(self.meta["area"])    # u.m2
        canvas_cube_hdu.data *= area.to(u.cm**2)

        bin_widths = np.diff(fov_waveset)
        bin_widths = 0.5 * (np.r_[0, bin_widths] + np.r_[bin_widths, 0])
        canvas_cube_hdu.data *= bin_widths.to(u.AA)

        return canvas_cube_hdu

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

    @property
    def waveset(self):
        """Returns a wavelength vector in um"""
        if len(self.cube_fields) > 0:
            i = np.argmax([cube.header["NAXIS3"] for cube in field_cubes])
            _waveset = fu.get_cube_waveset(field_cubes[i].header,
                                              return_quantity=True)
        else:
            i = np.argmax([len(spec.waveset) for spec in self.spectra.values()])
            _waveset = self.spectra[i].waveset

        return _waveset.to(u.um)

    @property
    def cube_fields(self):
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header["NAXIS"] == 3]

    @property
    def image_fields(self):
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header["NAXIS"] == 2]

    @property
    def table_fields(self):
        return [field for field in self.fields
                if isinstance(field, Table)]

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
