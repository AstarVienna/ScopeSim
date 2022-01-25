from copy import deepcopy

import numpy as np
from scipy.interpolate import interp1d

from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
from synphot import Empirical1D, SourceSpectrum
from synphot.units import PHOTLAM

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

    def __init__(self, header, waverange, detector_header=None, **kwargs):
        self.meta = {"id": None,
                     "wave_min": utils.quantify(waverange[0], u.um),
                     "wave_max": utils.quantify(waverange[1], u.um),
                     "wave_bin_n": 1,
                     "wave_bin_type": "linear",

                     "area": 0 * u.m**2,
                     "pixel_scale": "!INST.pixel_scale",
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

        self.header = fits.Header()
        self.header["NAXIS"] = 2
        self.header["NAXIS1"] = header["NAXIS1"]
        self.header["NAXIS2"] = header["NAXIS2"]
        self.header.update(header)
        self.detector_header = detector_header
        self.hdu = None

        self.image_plane_id = 0
        self.fields = []
        self.spectra = {}

        self.cube = None        # 3D array for IFU, long-lit, Slicer-MOS
        self.image = None       # 2D array for Imagers
        self.spectrum = None    # SourceSpectrum for Fibre-fed MOS

        self._waverange = None
        self._wavelength = None
        self._volume = None

    def extract_from(self, src):
        """ ..assumption: Bandpass has been applied

        .. note:: Spectra are cut and copied from the original Source object.
            They are in original units. ph/s/pix comes in the make_**** methods

        """

        if not isinstance(src, SourceBase):
            raise ValueError("source must be a Source object: {}"
                             "".format(type(src)))

        fields_in_fov = [field for field in src.fields
                         if fu.is_field_in_fov(self.header, field)]

        spec_refs = []
        volume = self.volume()
        for i in range(len(fields_in_fov)):
            fld = fields_in_fov[i]
            if isinstance(fld, Table):
                fields_in_fov[i] = fu.extract_area_from_table(fld, volume)
                spec_refs += list(np.unique(fields_in_fov[i] ["ref"]))

            elif isinstance(fld, fits.ImageHDU):
                if fld.header["NAXIS"] in (2, 3):
                    fields_in_fov[i] = fu.extract_area_from_imagehdu(fld, volume)
                if fld.header["NAXIS"] == 2 or fld.header.get("BG_SRC"):
                    ref = fld.header.get("SPEC_REF")
                    if ref is not None:
                        spec_refs += [ref]

        waves = volume["waves"] * u.Unit(volume["wave_unit"])
        spectra = {ref: fu.extract_range_from_spectrum(src.spectra[ref], waves)
                   for ref in np.unique(spec_refs)}

        self.fields = fields_in_fov
        self.spectra = spectra

    def view(self, hdu_type="image", sub_pixel=None, use_photlam=None):
        """
        Forces the self.fields to be viewed as a single object

        Parameters
        ----------
        sub_pixel : bool
        hdu_type : str
            ["cube", "image", "spectrum"]

        Returns
        -------
        self.hdu : fits.ImageHDU, synphot.SourceSpectrum

        """
        if sub_pixel is not None:
            self.meta["sub_pixel"] = sub_pixel

        if hdu_type == "image":
            use_photlam = False if use_photlam is None else use_photlam
            self.hdu = self.make_image_hdu(use_photlam=use_photlam)
        elif hdu_type == "cube":
            self.hdu = self.make_cube_hdu()
        elif hdu_type == "spectrum":
            self.hdu = self.make_spectrum()

        return self.hdu

    def flatten(self):
        if self.hdu.header["NAXIS"] == 3:
            image = np.sum(self.hdu.data, axis=0)
            self.hdu.data = image

    def make_spectrum(self):
        """
        This is needed for when we do incoherent MOS instruments.
        Each fibre doesn't care about the spatial information.

        1. Make waveset and zeros flux
            make dict of spectra evaluated at waveset

        2. Find Cube fields
            collapse cube along spatial dimensions --> spectrum vector
            convert vector to PHOTLAM
            interpolate at waveset
            add to zeros flux vector

        3. Find Image fields
            sum image over both dimensions
            evaluate SPEC_REF spectum at waveset
            multiply by sum
            add to zeros flux vector

        4. Find Table fields
            evaluate all spectra at waveset
            for each unique ref, sum the weights
            add each spectra * sum of weights to zeros flux vector

        Returns
        -------
        spec : SourceSpectrum
            [PHOTLAM]

        """
        fov_waveset = self.waveset
        canvas_flux = np.zeros(len(fov_waveset))

        specs = {ref: spec(fov_waveset) for ref, spec in self.spectra.items()}

        for field in self.cube_fields:
            hdu_waveset = fu.get_cube_waveset(field.header, return_quantity=True)
            fluxes = field.data.sum(axis=2).sum(axis=1)
            fov_waveset_fluxes = np.interp(fov_waveset, hdu_waveset, fluxes)

            field_unit = field.header.get("BUNIT", PHOTLAM)
            flux_scale_factor = u.Unit(field_unit).to(PHOTLAM)

            canvas_flux += fov_waveset_fluxes * flux_scale_factor

        for field in self.image_fields:
            ref = field.header["SPEC_REF"]
            weight = np.sum(field.data)
            canvas_flux += self.spectra[ref](fov_waveset).value * weight


        for field in self.table_fields:
            refs = np.array(field["ref"])
            weights = np.array(field["weight"])
            weight_sums = {ref: np.sum(weights[refs == ref])
                           for ref in np.unique(refs)}
            for ref, weight in weight_sums.items():
                canvas_flux += self.spectra[ref](fov_waveset).value * weight

        for field in self.background_fields:
            bg_solid_angle = u.Unit(field.header["SOLIDANG"]).to(u.arcsec**-2)
            pixel_area = utils.from_currsys(self.meta["pixel_scale"]) ** 2
            area_factor = pixel_area * bg_solid_angle       # arcsec**2 * arcsec**-2

            ref = field.header["SPEC_REF"]
            canvas_flux += self.spectra[ref](fov_waveset).value * area_factor

        spectrum = SourceSpectrum(Empirical1D, points=fov_waveset,
                                  lookup_table=canvas_flux)

        return spectrum

    def make_image_hdu(self, use_photlam=False):
        """
        Used for imaging

        Output image units are ph s-1 pixel-1

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

        Parameters
        ----------
        use_photlam : bool
            Default False. Defines the flux units of the image pixels

        Returns
        -------
        image_hdu : fits.ImageHDU
            [ph s-1 pixel-1] or PHOTLAM (if use_photlam=True)

        """
        spline_order = utils.from_currsys("!SIM.computing.spline_order")

        # 1. Make waveset and canvas image
        fov_waveset = self.waveset
        bin_widths = np.diff(fov_waveset)       # u.um
        bin_widths = 0.5 * (np.r_[0, bin_widths] + np.r_[bin_widths, 0])
        area = utils.from_currsys(self.meta["area"])    # u.m2

        # PHOTLAM * u.um * u.m2 --> ph / s
        specs = {ref: spec(fov_waveset) for ref, spec in self.spectra.items()}
        if use_photlam is False:
            for key in specs:
                specs[key] = (specs[key] * bin_widths * area).to(u.ph / u.s)

        fluxes = {ref: np.sum(spec).value for ref, spec in specs.items()}
        naxis1, naxis2 = self.header["NAXIS1"], self.header["NAXIS2"]
        canvas_image_hdu = fits.ImageHDU(data=np.zeros((naxis2, naxis1)),
                                         header=self.header)

        # 2. Find Cube fields
        for field in self.cube_fields:
            image = np.sum(field.data, axis=0)
            tmp_hdu = fits.ImageHDU(data=image, header=field.header)
            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(
                tmp_hdu,
                canvas_image_hdu,
                conserve_flux=True,
                spline_order=spline_order)

        # 2. Find Image fields
        for field in self.image_fields:
            image = deepcopy(field.data)
            tmp_hdu = fits.ImageHDU(data=image, header=field.header)
            tmp_hdu.data *= fluxes[field.header["SPEC_REF"]]  # ph / s
            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(
                tmp_hdu,
                canvas_image_hdu,
                conserve_flux=True,
                spline_order=spline_order)

        # 3. Find Table fields
        for field in self.table_fields:
            # x, y are ALWAYS in arcsec - crval is in deg
            xpix, ypix = imp_utils.val2pix(self.header,
                                           field["x"] / 3600,
                                           field["y"] / 3600)
            if utils.from_currsys(self.meta["sub_pixel"]):
                for i, row in enumerate(field):
                    xs, ys, fracs = imp_utils.sub_pixel_fractions(xpix[i],
                                                                  ypix[i])
                    ref, weight = row["ref"], row["weight"]
                    for x, y, f in zip(xs, ys, fracs):
                        canvas_image_hdu.data[y, x] += fluxes[ref] * weight * f
            else:
                x = np.array(xpix + 0.5).astype(int)
                y = np.array(ypix + 0.5).astype(int)     # quickest way to round
                f = np.array([fluxes[ref] for ref in field["ref"]])
                weight = np.array(field["weight"])
                canvas_image_hdu.data[y, x] += f * weight

        # 4. Find Background fields
        for field in self.background_fields:
            bg_solid_angle = u.Unit(field.header["SOLIDANG"]).to(u.arcsec**-2)
            pixel_area = utils.from_currsys(self.meta["pixel_scale"]) ** 2
            area_factor = pixel_area * bg_solid_angle       # arcsec**2 * arcsec**-2

            flux = fluxes[field.header["SPEC_REF"]] * area_factor
            canvas_image_hdu.data += flux

        image_hdu = canvas_image_hdu        # [ph s-1 pixel-1]

        return image_hdu

    def make_cube_hdu(self):
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
        - tables: (PHOTLAM in spectrum)
        - images: arcsec-2 (PHOTLAM in spectrum)
        - cubes: PHOTLAM arcsec-2

        .. warning:: Input Images and Cubes should have units of PHOTLAM arcsec-2

        Returns
        -------
        canvas_cube_hdu : fits.ImageHDU
            [ph s-1 AA-1 arcsec-2]      # as needed by SpectralTrace

        """
        spline_order = utils.from_currsys("!SIM.computing.spline_order")


        # 1. Make waveset and canvas cube (area, bin_width are applied at end)

        wave_min = self.meta["wave_min"]        # Quantity [um]
        wave_max = self.meta["wave_max"]

        wave_unit = u.Unit(utils.from_currsys("!SIM.spectral.wave_unit"))
        dwave = utils.from_currsys("!SIM.spectral.spectral_bin_width")  # Not a quantity
        fov_waveset = np.arange(wave_min.value, wave_max.value, dwave) * wave_unit
        fov_waveset = fov_waveset.to(u.um)

        # fov_waveset = self.waveset
        # wave_bin_n = len(fov_waveset)
        # if "lin" in self.meta["wave_bin_type"]:
        #     fov_waveset = np.linspace(wave_min, wave_max, wave_bin_n)
        # elif "log" in self.meta["wave_bin_type"]:
        #     wmin, wmax = wave_min.to(u.um).value, wave_max.to(u.um).value
        #     fov_waveset = np.logspace(wmin, wmax, wave_bin_n)

        specs = {ref: spec(fov_waveset)                     # PHOTLAM = ph/s/cm2/AA
                 for ref, spec in self.spectra.items()}

        # make canvas cube based on waveset of largest cube and NAXIS1,2 from fov.header
        naxis1, naxis2 = self.header["NAXIS1"], self.header["NAXIS2"]
        naxis3 = len(fov_waveset)
        canvas_cube_hdu = fits.ImageHDU(data=np.zeros((naxis3, naxis2, naxis1)),
                                        header=self.header)
        canvas_cube_hdu.header["BUNIT"] = "ph s-1 cm-2 AA-1"

        # 2. Add Cube fields
        for field in self.cube_fields:
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Assumption is that ImageHDUs have units of PHOTLAM arcsec-2
            # ..todo: Add a catch to get ImageHDU with BUNITs
            field_waveset = fu.get_cube_waveset(field.header,
                                                return_quantity=True)
            # ..todo: Deal with this bounds_error in a more elegant way
            field_interp = interp1d(field_waveset.to(u.um).value,
                                    field.data, axis=0, kind="linear",
                                    bounds_error=False, fill_value=0)
            field_data = field_interp(fov_waveset.value)
            field_unit = field.header.get("BUNIT", "ph s-1 cm-2 AA-1")
            eq = u.spectral_density(fov_waveset)
            flux_scale_factor = u.Unit(field_unit).to("ph s-1 cm-2 AA-1",
                                                      equivalencies=eq)

            # OC [2021-12-14] add extra dimensions for layer-wise multiplication of the cube
            field_data *= flux_scale_factor[:, None, None]
            field_hdu = fits.ImageHDU(data=field_data, header=field.header)
            canvas_cube_hdu = imp_utils.add_imagehdu_to_imagehdu(field_hdu,
                                                    canvas_cube_hdu,
                                                    spline_order=spline_order)

        # 3. Find Image fields
        for field in self.image_fields:
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Assumption is that ImageHDUs have units of PHOTLAM arcsec-2
            # ..todo: Add a catch to get ImageHDU with BUNITs
            canvas_image_hdu = fits.ImageHDU(data=np.zeros((naxis2, naxis1)),
                                             header=self.header)
            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(field,
                                                    canvas_image_hdu,
                                                    spline_order=spline_order)
            spec = specs[field.header["SPEC_REF"]]
            field_cube = canvas_image_hdu.data[None, :, :] * spec[:, None, None]  # 2D * 1D -> 3D
            canvas_cube_hdu.data += field_cube.value

        # 4. Find Table fields
        for field in self.table_fields:
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Point sources are in PHOTLAM per pixel
            # Point sources need to be scaled up by inverse pixel_area
            pixel_area = utils.from_currsys(self.meta["pixel_scale"]) ** 2
            for row in field:
                xsky, ysky = row["x"], row["y"]
                ref, weight = row["ref"], row["weight"]
                # x, y are ALWAYS in arcsec - crval is in deg
                xpix, ypix = imp_utils.val2pix(self.header, xsky / 3600, ysky / 3600)
                if utils.from_currsys(self.meta["sub_pixel"]):
                    xs, ys, fracs = imp_utils.sub_pixel_fractions(xpix, ypix)
                    for i, j, k in zip(xs, ys, fracs):
                        flux_vector = specs[ref].value * weight * k / pixel_area
                        canvas_cube_hdu.data[:, j, i] +=  flux_vector
                else:
                    x, y = int(xpix), int(ypix)
                    flux_vector = specs[ref].value * weight / pixel_area
                    canvas_cube_hdu.data[:, y, x] += flux_vector

        # 5. Add Background fields
        for field in self.background_fields:
            # bg_solid_angle = u.Unit(field.header["SOLIDANG"]).to(u.arcsec**-2)  # float [arcsec-2]
            # pixel_area = utils.from_currsys(self.meta["pixel_scale"]) ** 2      # float [arcsec2]
            # area_factor = pixel_area * bg_solid_angle                           # float [arcsec2 * arcsec-2]

            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # spec = specs[field.header["SPEC_REF"]] * area_factor
            spec = specs[field.header["SPEC_REF"]]
            canvas_cube_hdu.data += spec[:, None, None].value

        # 6. Convert from PHOTLAM to ph/s/voxel
        #    PHOTLAM = ph/s/cm-2/AA
        #    area = m2, fov_waveset = um
        # SpectralTrace wants ph/s/um/arcsec2 --> get rid of m2, leave um
        area = utils.from_currsys(self.meta["area"])  # u.m2
        canvas_cube_hdu.data *= area.to(u.cm ** 2).value
        canvas_cube_hdu.data *= 1e4        # ph/s/AA/arcsec2 --> ph/s/um/arcsec2

        # bin_widths = np.diff(fov_waveset).to(u.AA).value
        # bin_widths = 0.5 * (np.r_[0, bin_widths] + np.r_[bin_widths, 0])
        # canvas_cube_hdu.data *= bin_widths[:, None, None]

        cdelt3 = np.diff(fov_waveset[:2])[0]
        canvas_cube_hdu.header.update({"CDELT3": cdelt3.to(u.um).value,
                                       "CRVAL3": fov_waveset[0].value,
                                       "CRPIX3": 1,
                                       "CUNIT3": "um",
                                       "CTYPE3": "WAVE"})
        # ..todo: Add the log wavelength keyword here, if log scale is needed

        cube_hdu = canvas_cube_hdu      # [ph s-1 AA-1 (arcsec-2)]

        return cube_hdu

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
    def data(self):
        if self.hdu is not None:
            data = self.hdu.data
        elif self.image is not None:
            data = self.image
        elif self.cube is not None:
            data = self.cube
        elif self.spectrum is not None:
            data = self.spectrum
        else:
            data = None

        return data

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
        return utils.quantify(self._wavelength, u.um)

    @property
    def waveset(self):
        """Returns a wavelength vector in um"""

        field_cubes = self.cube_fields
        if len(field_cubes) > 0:
            i = np.argmax([cube.header["NAXIS3"] for cube in field_cubes])
            _waveset = fu.get_cube_waveset(field_cubes[i].header,
                                           return_quantity=True)
        elif len(self.spectra) > 0:
            wavesets = [spec.waveset for spec in self.spectra.values()]
            _waveset = np.unique(np.concatenate(wavesets))
        else:
            _waveset = self.waverange * u.um

        return _waveset.to(u.um)

    @property
    def cube_fields(self):
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header["NAXIS"] == 3
                and field.header.get("BG_SRC", False) is False]

    @property
    def image_fields(self):
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header["NAXIS"] == 2
                and field.header.get("BG_SRC", False) is False]

    @property
    def table_fields(self):
        return [field for field in self.fields
                if isinstance(field, Table)]


    @property
    def background_fields(self):
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header.get("BG_SRC", False) is True]

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
