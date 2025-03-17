"""Defines FieldOfView class."""

from copy import deepcopy
from itertools import chain
from collections.abc import Iterable

import numpy as np
from scipy.interpolate import interp1d

from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from synphot import Empirical1D, SourceSpectrum
from synphot.units import PHOTLAM

from . import fov_utils as fu
from . import image_plane_utils as imp_utils

from ..base_classes import SourceBase, FieldOfViewBase
from ..utils import from_currsys, quantify, has_needed_keywords, get_logger


logger = get_logger(__name__)


class FieldOfView(FieldOfViewBase):
    """
    A FOV is spectro-spatial volume cut out of a Source object.

    Flux units after extracting the fields from the Source are in ph/s/pixel

    The initial header should contain an on-sky WCS description:
    - CDELT-, CUNIT-, NAXIS- : for pixel scale and size (assumed CUNIT in deg)
    - CRVAL-, CRPIX- : for positioning the final image
    - CTYPE- : is assumed to be "RA---TAN", "DEC--TAN"

    and an image-plane WCS description
    - CDELT-D, CUNIT-D, NAXISn : for pixel scale and size (assumed CUNIT in mm)
    - CRVAL-D, CRPIX-D : for positioning the final image
    - CTYPE-D : is assumed to be "LINEAR", "LINEAR"

    The wavelength range is given by waverange

    """

    def __init__(self, header, waverange, detector_header=None, cmds=None, **kwargs):
        self.meta = {
            "id": None,
            "wave_min": quantify(waverange[0], u.um),
            "wave_max": quantify(waverange[1], u.um),
            "wave_bin_n": 1,
            "wave_bin_type": "linear",

            "area": 0 * u.m**2,
            "pixel_area": None,  # [arcsec]
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

        self.cmds = cmds

        if not any((has_needed_keywords(header, s) for s in ["", "S"])):
            raise ValueError(
                f"Header must contain a valid sky-plane WCS: {dict(header)}")
        if not has_needed_keywords(header, "D"):
            raise ValueError(
                f"Header must contain a valid image-plane WCS: {dict(header)}")

        self.header = fits.Header()
        self.header["NAXIS"] = 2
        self.header["NAXIS1"] = header["NAXIS1"]
        self.header["NAXIS2"] = header["NAXIS2"]
        self.header.update(header)
        self._ensure_deg_header()
        self.detector_header = detector_header
        self.hdu = None

        self.image_plane_id = 0
        self.fields = []
        self.spectra = {}

        # These are apparently not supposed to be used?
        self.cube = None        # 3D array for IFU, long-lit, Slicer-MOS
        # self.image = None       # 2D array for Imagers
        # self.spectrum = None    # SourceSpectrum for Fibre-fed MOS

        self._waverange = None
        self._wavelength = None
        self._volume = None

    def _pixarea(self, hdr):
        return (hdr["CDELT1"] * u.Unit(hdr["CUNIT1"]) *
                hdr["CDELT2"] * u.Unit(hdr["CUNIT2"])).to(u.arcsec ** 2)

    @property
    def trace_id(self):
        """Return the name of the trace."""
        return self.meta["trace_id"]

    @property
    def pixel_area(self):
        if self.meta["pixel_area"] is None:
            # [arcsec] (really?)
            self.meta["pixel_area"] = self._pixarea(self.header).value
        return self.meta["pixel_area"]

    def extract_from(self, src):
        """..assumption: Bandpass has been applied.

        .. note:: Spectra are cut and copied from the original Source object.
            They are in original units. ph/s/pix comes in the make_**** methods

        """
        if not isinstance(src, SourceBase):
            raise ValueError(f"source must be a Source object: {type(src)}")

        fields_in_fov = [field.field for field in src.fields
                         if fu.is_field_in_fov(self.header, field)]
        if not fields_in_fov:
            logger.warning("No fields in FOV.")

        spec_refs = set()
        volume = self.volume()
        for ifld, fld in enumerate(fields_in_fov):
            if isinstance(fld, Table):
                extracted_field = fu.extract_area_from_table(fld, volume)
                spec_refs.update(extracted_field["ref"])
                fields_in_fov[ifld] = extracted_field

            elif isinstance(fld, fits.ImageHDU):
                if fld.header["NAXIS"] in (2, 3):
                    extracted = fu.extract_area_from_imagehdu(fld, volume)
                    try:
                        fill_value = self.cmds["!SIM.computing.nan_fill_value"]
                    except TypeError:
                        # This occurs in testing, not sure why
                        fill_value = 0.
                    if np.ma.fix_invalid(extracted.data, copy=False,
                                         fill_value=fill_value).mask.any():
                        logger.warning("Source contained NaN values, which "
                                       "were replaced by %f.", fill_value)
                        logger.info(
                            "The replacement value for NaNs in sources can be "
                            "set in '!SIM.computing.nan_fill_value'.")
                    fields_in_fov[ifld] = extracted
                if fld.header["NAXIS"] == 2 or fld.header.get("BG_SRC"):
                    ref = fld.header.get("SPEC_REF")
                    if ref is not None:
                        spec_refs.add(ref)

        waves = volume["waves"] * u.Unit(volume["wave_unit"])
        spectra = {ref: fu.extract_range_from_spectrum(src.spectra[int(ref)], waves)
                   for ref in spec_refs}

        self.fields = fields_in_fov
        self.spectra = spectra

    def view(self, hdu_type="image", sub_pixel=None, use_photlam=None):
        """
        Force the self.fields to be viewed as a single object.

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
            # FIXME: why not just make False the default value??
            use_photlam = False if use_photlam is None else use_photlam
            self.hdu = self.make_image_hdu(use_photlam=use_photlam)
        elif hdu_type == "cube":
            self.hdu = self.make_cube_hdu()
        elif hdu_type == "spectrum":
            self.hdu = self.make_spectrum()

        return self.hdu

    def flatten(self):
        """If cube, collapse along first axis."""
        if self.hdu and self.hdu.header["NAXIS"] == 3:
            image = np.sum(self.hdu.data, axis=0)
            self.hdu.data = image

    def _evaluate_spectrum_with_weight(self, ref, waveset, weight):
        return self.spectra[int(ref)](waveset).value * weight

    def _calc_area_factor(self, field):
        bg_solid_angle = u.Unit(field.header["SOLIDANG"]).to(u.arcsec**-2)
        # arcsec**2 * arcsec**-2
        area_factor = self.pixel_area * bg_solid_angle
        return area_factor

    def _make_spectrum_cubefields(self, fov_waveset):
        """
        Find Cube fields.

        * collapse cube along spatial dimensions --> spectrum vector
        * convert vector to PHOTLAM
        * interpolate at waveset
        * yield scaled flux to be added to canvas flux
        """
        for field in self.cube_fields:
            hdu_waveset = fu.get_cube_waveset(field.header,
                                              return_quantity=True)
            fluxes = field.data.sum(axis=2).sum(axis=1)
            fov_waveset_fluxes = np.interp(fov_waveset, hdu_waveset, fluxes)

            field_unit = field.header.get("BUNIT", PHOTLAM)
            flux_scale_factor = u.Unit(field_unit).to(PHOTLAM)

            yield fov_waveset_fluxes * flux_scale_factor

    def _make_spectrum_imagefields(self, waveset):
        """
        Find Image fields.

        * sum image over both dimensions
        * evaluate SPEC_REF spectum at waveset
        * yield spectrum multiply by sum to be added to canvas flux
        """
        for field in self.image_fields:
            weight = np.sum(field.data)
            yield self._evaluate_spectrum_with_weight(field.header["SPEC_REF"],
                                                      waveset, weight)

    def _make_spectrum_tablefields(self, waveset):
        """
        Find Table fields.

        * evaluate all spectra at waveset
        * for each unique ref, sum the weights
        * yield each spectrum * sum of weights to be added to canvas flux
        """
        for field in self.table_fields:
            refs = np.array(field["ref"])
            weights = np.array(field["weight"])
            # TODO: could do grouping of table with both columns??
            for ref in set(refs):
                weight = np.sum(weights, where=refs == ref)
                yield self._evaluate_spectrum_with_weight(ref, waveset, weight)

    def _make_spectrum_backfields(self, waveset):
        for field in self.background_fields:
            yield self._evaluate_spectrum_with_weight(
                field.header["SPEC_REF"], waveset,
                self._calc_area_factor(field))

    def make_spectrum(self):
        """
        TBA.

        This is needed for when we do incoherent MOS instruments.
        Each fibre doesn't care about the spatial information.

        Returns
        -------
        spec : SourceSpectrum
            [PHOTLAM]

        """
        fov_waveset = self.waveset
        # Start with zero flux no ensure correct array shape even if none of
        # the sub-functions yield anything.
        canvas_flux = sum(chain(
            self._make_spectrum_cubefields(fov_waveset),
            self._make_spectrum_imagefields(fov_waveset),
            self._make_spectrum_tablefields(fov_waveset),
            self._make_spectrum_backfields(fov_waveset),
        ), start=np.zeros_like(fov_waveset.value))

        spectrum = SourceSpectrum(Empirical1D, points=fov_waveset,
                                  lookup_table=canvas_flux)
        return spectrum

    def _make_image_cubefields(self, area):
        """
        Find Cube fields.

        * collapse cube along wavelength axis
        * rescale and reorient image
        * yield cube image  to be added to canvas image
        """
        for field in self.cube_fields:
            # cube_fields come in with units of photlam/arcsec2,
            # need to convert to ph/s
            # We need to the voxel volume (spectral and solid angle) for that.
            # ..todo: implement branch for use_photlam is True
            spectral_bin_width = (field.header["CDELT3"] *
                                  u.Unit(field.header["CUNIT3"])
                                  ).to(u.Angstrom)
            # First collapse to image, then convert units
            image = np.sum(field.data, axis=0) * PHOTLAM/u.arcsec**2
            image = (image * self._pixarea(field.header) * area *
                     spectral_bin_width).to(u.ph/u.s)
            yield fits.ImageHDU(data=image, header=field.header)

    def _make_image_imagefields(self, fluxes):
        """
        Find Image fields.

        * sum spectra between wavelength edges
        * multiply image by summed spectrum
        * yield image  to be added to canvas image
        """
        for field in self.image_fields:
            image = deepcopy(field.data)
            image *= fluxes[field.header["SPEC_REF"]]  # ph / s
            yield fits.ImageHDU(data=image, header=field.header)

    def _make_image_tablefields(self, fluxes):
        """
        Find Table fields.

        * sum spectra between wavelength edges
        * yield summed flux at x,y position to be added to canvas image
        """
        for field in self.table_fields:
            # x, y are ALWAYS in arcsec - crval is in deg
            xpix, ypix = imp_utils.val2pix(self.header,
                                           field["x"] / 3600,
                                           field["y"] / 3600)
            if from_currsys(self.meta["sub_pixel"], self.cmds):
                for idx, row in enumerate(field):
                    xs, ys, fracs = imp_utils.sub_pixel_fractions(xpix[idx],
                                                                  ypix[idx])
                    for x, y, frac in zip(xs, ys, fracs):
                        yield fluxes[row["ref"]] * frac, row["weight"], x, y
            else:
                # Note: these had x/ypix+0.5 until a06ab75
                # TODO: could these be something more numpythonic grid-ish?
                x = np.array(xpix).astype(int)
                y = np.array(ypix).astype(int)     # quickest way to round
                flux = np.array([fluxes[int(ref)] for ref in field["ref"]])
                yield flux, np.array(field["weight"]), x, y

    def _make_image_backfields(self, fluxes):
        for field in self.background_fields:
            flux = fluxes[field.header["SPEC_REF"]]
            yield flux * self._calc_area_factor(field)

    def make_image_hdu(self, use_photlam=False):
        """
        TBA.

        Used for imaging.

        Output image units are ph s-1 pixel-1

        .. note:: ``self.make_image()`` does NOT store anything in ``self.image``

            See make_cube for an explanation

        Make canvas image from NAXIS1,2 from fov.header

        Parameters
        ----------
        use_photlam : bool
            Default False. Defines the flux units of the image pixels

        Returns
        -------
        image_hdu : fits.ImageHDU
            [ph s-1 pixel-1] or PHOTLAM (if use_photlam=True)

        """
        spline_order = from_currsys("!SIM.computing.spline_order", self.cmds)
        # Make waveset and canvas image
        fov_waveset = self.waveset
        bin_widths = np.diff(fov_waveset)       # u.um
        bin_widths = 0.5 * (np.r_[0, bin_widths] + np.r_[bin_widths, 0])
        area = from_currsys(self.meta["area"], self.cmds) << u.m**2   # u.m2

        # PHOTLAM * u.um * u.m2 --> ph / s
        specs = {ref: spec(fov_waveset) if use_photlam
                 else (spec(fov_waveset) * bin_widths * area).to(u.ph / u.s)
                 for ref, spec in self.spectra.items()}

        fluxes = {ref: np.sum(spec.value) for ref, spec in specs.items()}

        canvas_image_hdu = fits.ImageHDU(
            data=np.zeros((self.header["NAXIS2"], self.header["NAXIS1"])),
            header=self.header)

        for tmp_hdu in chain(self._make_image_cubefields(area),
                             self._make_image_imagefields(fluxes)):
            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(
                tmp_hdu,
                canvas_image_hdu,
                conserve_flux=True,
                spline_order=spline_order)

        for flux, weight, x, y in self._make_image_tablefields(fluxes):
            if from_currsys(self.meta["sub_pixel"], self.cmds):
                # These x and y should not be arrays when sub_pixel is
                # enabled, it is therefore not necessary to deploy the fix
                # below in the else-branch.
                assert not isinstance(x, Iterable), "x must be an integer"
                try:
                    canvas_image_hdu.data[y, x] += flux * weight
                except IndexError:
                    continue
            else:
                # Mask out any stars that were pushed out of the fov by rounding
                mask = ((x < canvas_image_hdu.data.shape[1]) *
                        (y < canvas_image_hdu.data.shape[0]))

                # This used to contain this line:
                #   canvas_image_hdu.data[y[mask], x[mask]] += flux[mask] * weight[mask]
                # However, that is wrong when there are duplicate (x, y) pairs.
                # In those cases, only the last source flux is added to the
                # pixel. Therefor it is necessary to iterate over the sources.
                # The stacking of stars is tested by TestStackedStars in
                # test_flux_is_conserved_through_full_system.py

                for yi, xi, fluxi, weighti in zip(
                        y[mask], x[mask], flux[mask], weight[mask]):
                    canvas_image_hdu.data[yi, xi] += fluxi * weighti


        canvas_image_hdu.data = sum(self._make_image_backfields(fluxes),
                                    start=canvas_image_hdu.data)
        return canvas_image_hdu  # [ph s-1 pixel-1]

    def _make_cube_cubefields(self, fov_waveset):
        """
        Find Cube fields.

        * rescale and reorient cubes
        * interp1d smaller cubes with waveset
        * yield cubes to be added to cavas cube
        """
        for field in self.cube_fields:
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Assumption is that ImageHDUs have units of PHOTLAM arcsec-2
            field_waveset = fu.get_cube_waveset(field.header,
                                                return_quantity=True)

            # ..todo: Deal with this bounds_error in a more elegant way
            field_interp = interp1d(field_waveset.to(u.um).value,
                                    field.data, axis=0, kind="linear",
                                    bounds_error=False, fill_value=0)

            field_data = field_interp(fov_waveset.value)

            # Pixel scale conversion
            field_data *= self._pixarea(field.header).value / self.pixel_area
            field_hdu = fits.ImageHDU(data=field_data, header=field.header)
            yield field_hdu

    def _make_cube_imagefields(self, specs, spline_order):
        """
        Find Image fields.

        * rescale and reorient images
        * evaluate spectra at waveset
        * expand image by spectra to 3D form
        * yield image cubes to be added to cavas cube
        """
        for field in self.image_fields:
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Assumption is that ImageHDUs have units of PHOTLAM arcsec-2
            # ImageHDUs have photons/second/pixel.
            # ..todo: Add a catch to get ImageHDU with BUNITs
            canvas_image_hdu = fits.ImageHDU(
                data=np.zeros((self.header["NAXIS2"], self.header["NAXIS1"])),
                header=self.header)
            # FIX: Do not scale source data - make a copy first.
            # FIX: Use "Pixel scale conversion" as above.
            field_data = field.data * self._pixarea(field.header).value / self.pixel_area
            field_hdu = fits.ImageHDU(data=field_data, header=field.header)

            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(
                field_hdu,
                canvas_image_hdu,
                spline_order=spline_order)
            spec = specs[field.header["SPEC_REF"]]

            # 2D * 1D -> 3D
            field_cube = canvas_image_hdu.data[None, :, :] * spec[:, None, None]
            yield field_cube.value

    def _make_cube_tablefields(self, specs):
        """
        Find Table fields.

        * evaluate spectra at waveset
        * yield spectrum at x,y position to be added cavas to cube
        """
        for field in self.table_fields:
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Point sources are in PHOTLAM per pixel
            # Point sources need to be scaled up by inverse pixel_area
            for row in field:
                xsky, ysky = row["x"] / 3600, row["y"] / 3600
                # x, y are ALWAYS in arcsec - crval is in deg
                xpix, ypix = imp_utils.val2pix(self.header, xsky, ysky)
                flux_vector = specs[row["ref"]].value * row["weight"] / self.pixel_area

                if from_currsys(self.meta["sub_pixel"], self.cmds):
                    xs, ys, fracs = imp_utils.sub_pixel_fractions(xpix, ypix)
                    for i, j, k in zip(xs, ys, fracs):
                        yield flux_vector * k, i, j
                else:
                    yield flux_vector, int(xpix), int(ypix)

    def _make_cube_backfields(self, specs):
        for field in self.background_fields:
            # TODO: The following would have been identical to the other two
            #       make methods, but was commented out. Why?
            # bg_solid_angle = u.Unit(field.header["SOLIDANG"]).to(u.arcsec**-2)  # float [arcsec-2]
            # pixel_area = from_currsys(self.meta["pixel_scale"]) ** 2      # float [arcsec2]
            # area_factor = pixel_area * bg_solid_angle                           # float [arcsec2 * arcsec-2]

            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # spec = specs[field.header["SPEC_REF"]] * area_factor
            spec = specs[field.header["SPEC_REF"]]
            yield spec[:, None, None].value

    def make_cube_hdu(self):
        """
        TBA.

        Used for IFUs, slit spectrographs, and coherent MOSs (e.g.KMOS)

        Returned cube units are ``ph s-1 voxel-1``

        .. note:: ``self.make_cube()`` does NOT store anything in ``self.cube``

            self.cube and self.make_cube() are deliberately kept seperately
            so that self.cube will not be accidently overwritten by a rogue
            call from an Effect object.

            All Effect objects should specifically test whether
            ``self.cube is None`` before assigning a new cube it

        The cube is made with these steps:

        1. Make waveset and canvas cube::

            if at least one cube:
                set waveset to equal largest cube waveset
            else:
                make waveset from self.meta values
            make canvas cube based on waveset of largest cube and NAXIS1,2 from fov.header

        2. Find Cube fields (see ``FieldOfView._make_cube_cubefields()``).
        3. Find Image fields (see ``FieldOfView._make_cube_imagefields()``).
        4. Find Table fields (see ``FieldOfView._make_cube_tablefields()``).

        ``PHOTLAM = ph / (s * m2 * um)``.
        Original source fields are in units of:

        - tables: (PHOTLAM in spectrum)
        - images: arcsec-2 (PHOTLAM in spectrum)
        - cubes: PHOTLAM arcsec-2

        .. warning:: Input Images and Cubes should have units of PHOTLAM arcsec-2

        Returns
        -------
        canvas_cube_hdu : fits.ImageHDU
            [ph s-1 AA-1 arcsec-2]      # as needed by SpectralTrace

        """
        spline_order = from_currsys("!SIM.computing.spline_order", self.cmds)

        # 1. Make waveset and canvas cube (area, bin_width are applied at end)
        # TODO: Why is this not self.waveset? What's different?
        wave_unit = u.Unit(from_currsys("!SIM.spectral.wave_unit", self.cmds))
        fov_waveset = np.arange(
            self.meta["wave_min"].value, self.meta["wave_max"].value,
            from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)) * wave_unit
        fov_waveset = fov_waveset.to(u.um)

        # TODO: what's with this code??
        # fov_waveset = self.waveset
        # wave_bin_n = len(fov_waveset)
        # if "lin" in self.meta["wave_bin_type"]:
        #     fov_waveset = np.linspace(wave_min, wave_max, wave_bin_n)
        # elif "log" in self.meta["wave_bin_type"]:
        #     wmin, wmax = wave_min.to(u.um).value, wave_max.to(u.um).value
        #     fov_waveset = np.logspace(wmin, wmax, wave_bin_n)

        specs = {ref: spec(fov_waveset)                 # PHOTLAM = ph/s/cm2/AA
                 for ref, spec in self.spectra.items()}

        # make canvas cube based on waveset of largest cube and NAXIS1,2
        # from fov.header
        canvas_cube_hdu = fits.ImageHDU(
            data=np.zeros((len(fov_waveset),
                           self.header["NAXIS2"],
                           self.header["NAXIS1"])),
            header=self.header)
        canvas_cube_hdu.header["BUNIT"] = "ph s-1 cm-2 AA-1"

        for field_hdu in self._make_cube_cubefields(fov_waveset):
            canvas_cube_hdu = imp_utils.add_imagehdu_to_imagehdu(
                field_hdu,
                canvas_cube_hdu,
                spline_order=spline_order)

        canvas_cube_hdu.data = sum(self._make_cube_imagefields(specs,
                                                               spline_order),
                                   start=canvas_cube_hdu.data)

        for flux, x, y in self._make_cube_tablefields(specs):
            # To prevent adding array values in this manner.
            assert not isinstance(x, Iterable), "x should be integer"
            canvas_cube_hdu.data[:, y, x] += flux

        canvas_cube_hdu.data = sum(self._make_cube_backfields(specs),
                                   start=canvas_cube_hdu.data)

        # TODO: what's with this code??
        # 6. Convert from PHOTLAM to ph/s/voxel
        #    PHOTLAM = ph/s/cm-2/AA
        #    area = m2, fov_waveset = um
        # SpectralTrace wants ph/s/um/arcsec2 --> get rid of m2, leave um
        area = from_currsys(self.meta["area"], self.cmds)  # u.m2
        canvas_cube_hdu.data *= area.to(u.cm ** 2).value
        canvas_cube_hdu.data *= 1e4       # ph/s/AA/arcsec2 --> ph/s/um/arcsec2

        # TODO: what's with this code??
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
        return canvas_cube_hdu      # [ph s-1 AA-1 (arcsec-2)]

    def volume(self, wcs_prefix=""):
        xy = imp_utils.calc_footprint(self.header, wcs_suffix=wcs_prefix)
        unit = self.header[f"CUNIT1{wcs_prefix}"].lower()
        # FIXME: This is unused!!
        # wave_corners = self.waverange
        minmax = np.array((xy.min(axis=0), xy.max(axis=0)))
        self._volume = {"xs": minmax[:, 0],
                        "ys": minmax[:, 1],
                        "waves": self.waverange,
                        "xy_unit": unit,
                        "wave_unit": "um"}
        return self._volume

    @property
    def data(self):
        """Return either hdu.data, image, cube, spectrum or None."""
        if self.hdu is not None:
            return self.hdu.data
        if self.image is not None:
            return self.image
        if self.cube is not None:
            return self.cube
        if self.spectrum is not None:
            return self.spectrum
        return None

    @property
    def corners(self):
        """Return sky footprint, image plane footprint."""
        # Couldn't find where this is ever used, put warning here just in case
        logger.warning("calc_footprint has been updated, any code that "
                       "relies on this .corners property must likely be "
                       "adapted as well.")
        sky_corners = imp_utils.calc_footprint(self.header)
        imp_corners = imp_utils.calc_footprint(self.header, "D")
        return sky_corners, imp_corners

    @property
    def waverange(self):
        """Return wavelength range in um [wave_min, wave_max]."""
        if self._waverange is None:
            wave_min = quantify(self.meta["wave_min"], u.um).value
            wave_max = quantify(self.meta["wave_max"], u.um).value
            self._waverange = [wave_min, wave_max]
        return self._waverange

    @property
    def wavelength(self):
        """Return central wavelength in um."""
        if self._wavelength is None:
            self._wavelength = np.average(self.waverange)
        return quantify(self._wavelength, u.um)

    @property
    def waveset(self):
        """Return a wavelength vector in um."""
        if field_cubes := self.cube_fields:
            naxis3_max = np.argmax([cube.header["NAXIS3"]
                                    for cube in field_cubes])
            _waveset = fu.get_cube_waveset(field_cubes[naxis3_max].header,
                                           return_quantity=True)
        elif self.spectra:
            wavesets = [spec.waveset for spec in self.spectra.values()]
            _waveset = np.concatenate(wavesets)
        else:
            _waveset = self.waverange * u.um

        # TODO: tie the round to a global precision setting (this is defined
        #       better in another TODO somewhere...)
        # The rounding is necessary to not end up with things like:
        #   0.7 um
        #   0.7000000000000001 um
        #   0.7000000000000002 um
        # and yes, that actually happend...
        _waveset = np.unique(_waveset.to(u.um).round(10))
        return _waveset

    @property
    def cube_fields(self):
        """Return list of non-BG_SRC ImageHDU fields with NAXIS=3."""
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header["NAXIS"] == 3
                and not field.header.get("BG_SRC", False)]

    @property
    def image_fields(self):
        """Return list of non-BG_SRC ImageHDU fields with NAXIS=2."""
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header["NAXIS"] == 2
                and not field.header.get("BG_SRC", False)]

    @property
    def table_fields(self):
        """Return list of Table fields."""
        return [field for field in self.fields
                if isinstance(field, Table)]

    @property
    def background_fields(self):
        """Return list of BG_SRC ImageHDU fields."""
        return [field for field in self.fields
                if isinstance(field, fits.ImageHDU)
                and field.header.get("BG_SRC", False)]

    def _ensure_deg_header(self):
        cunit = u.Unit(self.header["CUNIT1"].lower())
        convf = cunit.to(u.deg)
        self.header["CDELT1"] *= convf
        self.header["CDELT2"] *= convf
        self.header["CRVAL1"] *= convf
        self.header["CRVAL2"] *= convf
        self.header["CUNIT1"] = "deg"
        self.header["CUNIT2"] = "deg"

    def __repr__(self):
        waverange = [self.meta["wave_min"].value, self.meta["wave_max"].value]
        msg = (f"{self.__class__.__name__}({self.header!r}, {waverange!r}, "
               f"{self.detector_header!r}, **{self.meta!r})")
        return msg

    def __str__(self):
        msg = (f"FOV id: {self.meta['id']}, with dimensions "
               f"({self.header['NAXIS1']}, {self.header['NAXIS2']})\n"
               f"Sky centre: ({self.header['CRVAL1']}, "
               f"{self.header['CRVAL2']})\n"
               f"Image centre: ({self.header['CRVAL1D']}, "
               f"{self.header['CRVAL2D']})\n"
               f"Wavelength range: ({self.meta['wave_min']}, "
               f"{self.meta['wave_max']})um\n")
        return msg
