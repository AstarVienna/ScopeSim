# -*- coding: utf-8 -*-
"""Defines FieldOfView class."""

from warnings import warn
from copy import deepcopy
from itertools import chain
from collections.abc import Iterable, Generator

import numpy as np
from scipy.interpolate import interp1d

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from synphot import Empirical1D, SourceSpectrum
from synphot.units import PHOTLAM

from . import image_plane_utils as imp_utils
from ..source.source_fields import (
    SourceField,
    SpectrumSourceField,
    ImageSourceField,
    CubeSourceField,
    TableSourceField,
    BackgroundSourceField,
)

from ..utils import (
    from_currsys,
    quantify,
    has_needed_keywords,
    get_logger,
    array_minmax,
    close_loop,
    unit_includes_per_physical_type,
)
from ..source.source import Source


logger = get_logger(__name__)


class FieldOfView:
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

    The wavelength range is given by waverange.

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
            "distortion": {
                "scale": [1, 1],
                "offset": [0, 0],
                "shear": [1, 1],
                "rotation": 0,
                "radius_of_curvature": None,
            },
            "conserve_image": True,
            "trace_id": None,
            "aperture_id": None,
        }
        self.meta.update(kwargs)

        self.cmds = cmds

        if not any(has_needed_keywords(header, s) for s in {"", "S"}):
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

        self.fields: list[SourceField] = []

        # These are apparently not supposed to be used?
        self.cube = None        # 3D array for IFU, long-lit, Slicer-MOS
        # self.image = None       # 2D array for Imagers
        # self.spectrum = None    # SourceSpectrum for Fibre-fed MOS

        self._waverange = None
        self._wavelength = None
        self._volume = None

    @staticmethod
    def _pixarea(hdr):
        return (hdr["CDELT1"] * u.Unit(hdr["CUNIT1"]) *
                hdr["CDELT2"] * u.Unit(hdr["CUNIT2"])).to(u.arcsec ** 2)

    @property
    def trace_id(self):
        """Return the name of the trace."""
        return self.meta["trace_id"]

    @property
    def pixel_area(self):
        """Return the area in arcsec**2 covered by one pixel."""
        if self.meta["pixel_area"] is None:
            # [arcsec] (really?)
            self.meta["pixel_area"] = self._pixarea(self.header).value
        return self.meta["pixel_area"]

    def extract_from(self, src) -> None:
        """
        Extract relevent fields from source object.

        Parameters
        ----------
        src : Source
            Input Source object to be "observed".

        Returns
        -------
        None

        Notes
        -----
        Spectra are cut and copied from the original Source object.
        They are in original units. ph/s/pix comes in the make_**** methods.

        This method assumes that Bandpass has been applied.

        """
        assert isinstance(src, Source), f"expected Source: {type(src)}"

        fields_in_fov = list(self.get_fields_in_fov(src.fields))

        if not fields_in_fov:
            logger.warning("No fields in FOV.")
        else:
            logger.debug("%d fields in FOV", len(fields_in_fov))

        corners_arcsec, _ = self.get_corners("arcsec")
        corners_deg, _ = self.get_corners("deg")
        minmax = array_minmax(corners_arcsec) * u.arcsec

        for field in fields_in_fov:
            if isinstance(field, TableSourceField):
                extracted = self.extract_area_from_table(field.field, minmax)
                if not len(extracted):
                    continue  # Can happen for small FOV and sparse table
                # TODO: Rework extract_area_from_table to also affect spectra
                #       and just return new copy of field.
                new_fld = TableSourceField(
                    field=extracted,
                    spectra={
                        ref: extract_range_from_spectrum(spec, self.waverange)
                        for ref, spec in field.spectra.items()
                        if ref in extracted["ref"]
                    },
                )
                self.fields.append(new_fld)

            elif isinstance(field, ImageSourceField):
                assert field.header["NAXIS"] == 2, "Invalid image HDU"
                extracted = self.extract_area_from_imagehdu(field.field, corners_deg)
                replace_nans(extracted, self.cmds)
                new_fld = ImageSourceField(
                    field=extracted,
                    spectra={
                        ref: extract_range_from_spectrum(spec, self.waverange)
                        for ref, spec in field.spectra.items()
                    },
                )
                self.fields.append(new_fld)

            elif isinstance(field, CubeSourceField):
                assert field.header["NAXIS"] == 3, "Invalid cube HDU"
                extracted = self.extract_area_from_imagehdu(field.field, corners_deg)
                replace_nans(extracted, self.cmds)
                new_fld = CubeSourceField(field=extracted)
                self.fields.append(new_fld)

            elif isinstance(field, BackgroundSourceField):
                new_fld = BackgroundSourceField(
                    field=None,
                    header=field.header,
                    spectra={
                        ref: extract_range_from_spectrum(spec, self.waverange)
                        for ref, spec in field.spectra.items()
                    },
                )
                self.fields.append(new_fld)

            else:
                raise TypeError(f"Unexpected source field type {type(field)}.")

    def view(
        self,
        hdu_type: str = "image",
        sub_pixel: bool | None = None,
        use_photlam: bool | None = None,
    ):
        """
        Force the self.fields to be viewed as a single object.

        Parameters
        ----------
        hdu_type : {"image", "cube", "spectrum"}
            DESCRIPTION.
        sub_pixel : bool | None, optional
            If None (the default), use value from meta.
        use_photlam : bool | None, optional
            If None (the default), assume False. Only used in imaging (why?).

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

    def get_fields_in_fov(self, fields: Iterable[SourceField]) -> Generator:
        """Return True if Source.field footprint is inside FOV footprint."""
        fov_corners, _ = self.get_corners("arcsec")

        for field in fields:
            field_corners = field.get_corners("arcsec")
            is_inside_fov = (
                (field_corners.max(axis=0) > fov_corners.min(axis=0)).all() and
                (field_corners.min(axis=0) < fov_corners.max(axis=0)).all()
            )
            if is_inside_fov:
                yield field

    @staticmethod
    def extract_area_from_table(table, minmax):
        """
        Extract the entries of a ``Table`` that fit inside the FOV volume.

        Parameters
        ----------
        table : table.Table
            The field table.
        minmax : quantity
            From FOV corners in the form of [[xmin, ymin], [xmax, ymax]].

        Returns
        -------
        cut_table : table.Table
            Table reduced to sources inside the FOV.

        """
        mask = ((table["x"].quantity >= minmax[0, 0]) *
                (table["x"].quantity < minmax[1, 0]) *
                (table["y"].quantity >= minmax[0, 1]) *
                (table["y"].quantity < minmax[1, 1]))
        return table[mask]

    def extract_area_from_imagehdu(self, imagehdu, corners):
        """
        Extract the part of a ``ImageHDU`` that fits inside the FOV volume.

        Parameters
        ----------
        imagehdu : fits.ImageHDU
            The field ImageHDU, either an image or a cube with wavelength [um].
        corners : quantity
            From FOV corners in deg.

        Returns
        -------
        new_imagehdu : fits.ImageHDU

        """
        # TODO: At some point streamline the debug logging here...
        hdr = imagehdu.header
        image_wcs = WCS(hdr, naxis=2)
        naxis1, naxis2 = hdr["NAXIS1"], hdr["NAXIS2"]
        logger.debug("old naxis: %s", [naxis1, naxis2])
        xy_hdu = image_wcs.calc_footprint(center=False, axes=(naxis1, naxis2))

        if image_wcs.wcs.cunit[0] == "deg":
            logger.debug("Found 'deg' in image WCS, applying 360 fix.")
            imp_utils._fix_360(xy_hdu)
        elif image_wcs.wcs.cunit[0] == "arcsec":
            logger.debug("Found 'arcsec' in image WCS, converting to deg.")
            xy_hdu *= u.arcsec.to(u.deg)

        logger.debug("XY HDU:\n%s", xy_hdu)
        logger.debug("XY FOV:\n%s", corners)

        xy0s = np.array((xy_hdu.min(axis=0), corners.min(axis=0))).max(axis=0)
        xy1s = np.array((xy_hdu.max(axis=0), corners.max(axis=0))).min(axis=0)
        logger.debug("xy0s: %s; xy1s: %s", xy0s, xy1s)

        # Round to avoid floating point madness
        xyp = image_wcs.wcs_world2pix(np.array([xy0s, xy1s]), 0).round(7)

        # To deal with negative CDELTs
        logger.debug("xyp:\n%s", xyp)
        xyp.sort(axis=0)
        logger.debug("xyp:\n%s", xyp)

        xy0p = np.max(((0, 0), np.floor(xyp[0]).astype(int)), axis=0)
        xy1p = np.min(((naxis1, naxis2), np.ceil(xyp[1]).astype(int)), axis=0)
        logger.debug("xy0p: %s; xy1p: %s", xy0p, xy1p)

        # Add 1 if the same
        xy1p += (xy0p == xy1p)
        logger.debug("xy0p: %s; xy1p: %s", xy0p, xy1p)

        new_wcs, new_naxis = imp_utils.create_wcs_from_points(
            np.array([xy0s, xy1s]).round(11), pixel_scale=hdr["CDELT1"])

        # TODO: Come back at some point and figure out if the failing tests
        #       here are relevant or can be ignored...
        # FIXME: Commented out for now because it appears too often IRL...
        # try:
        #     roundtrip = new_wcs.wcs_world2pix(
        #         np.array([xy0s, xy1s - .5*image_wcs.wcs.cdelt]), 0).round(5)
        #     np.testing.assert_array_equal(roundtrip[0], [0, 0])
        #     np.testing.assert_array_equal(roundtrip[1], new_naxis - [1, 1])
        # except AssertionError:
        #     logger.exception("WCS roundtrip assertion failed.")
        # FIXME: Related to the above, this sometimes fails:
        # np.testing.assert_equal(xy1p - xy0p, new_naxis)
        # This occurs when the floor and ceil in xy0s and xy1s produce an
        # off-by-one error. Using .round instead for both would solve things
        # in those cases, but breaks in other cases. Find a proper solution!
        # Note: This is not super fatal, because the resulting projections
        #       will trim off that extra pixel later on, but this should still
        #       be addressed.

        new_hdr = new_wcs.to_header()
        new_hdr.update({"NAXIS1": new_naxis[0], "NAXIS2": new_naxis[1]})

        if hdr["NAXIS"] == 3:
            new_hdr, data = self._extract_volume_from_cube(imagehdu, new_hdr,
                                                           xy0p, xy1p)
        else:
            data = imagehdu.data[xy0p[1]:xy1p[1],
                                 xy0p[0]:xy1p[0]]
            # new_hdr["SPEC_REF"] = hdr.get("SPEC_REF")

        if not data.size:
            logger.warning("Empty image HDU.")

        return fits.ImageHDU(header=new_hdr, data=data)

    def _extract_volume_from_cube(self, cubehdu, new_hdr, xy0p, xy1p):
        hdr = cubehdu.header
        # Look 0.5*wdel past the fov edges in each direction to catch any
        # slices where the middle wavelength value doesn't fall inside the
        # fov waverange, but up to 50% of the slice is actually inside the
        # fov waverange:
        # E.g. FOV: [1.92, 2.095], HDU bin centres: [1.9, 2.0, 2.1]
        # CDELT3 = 0.1, and HDU bin edges: [1.85, 1.95, 2.05, 2.15]
        # So 1.9 slice needs to be multiplied by 0.3, and 2.1 slice should be
        # multipled by 0.45 to reach the scaled contribution of the edge slices
        # This scaling factor is:
        # f = ((hdu_bin_centre - fov_edge [+/-] 0.5 * cdelt3) % cdelt3) / cdelt3

        swcs = WCS(hdr).spectral
        with u.set_enabled_equivalencies(u.spectral()):
            hdu_waves = swcs.pixel_to_world(np.arange(swcs.pixel_shape[0]))
            hdu_waves = hdu_waves.quantity.to(u.um)

        wdel = hdr["CDELT3"] * u.Unit(hdr.get("CUNIT3", "AA"))
        fov_wmin, fov_wmax = self.waverange
        # need to go [+/-] half a bin
        mask = ((hdu_waves > fov_wmin - 0.5 * wdel) *
                (hdu_waves <= fov_wmax + 0.5 * wdel))

        # OC [2021-12-14] if fov range is not covered by the source return nothing
        if not np.any(mask):
            logger.warning("FOV %s um - %s um: not covered by Source",
                           fov_wmin, fov_wmax)
            # FIXME: returning None here breaks the principle that a function
            #        should always return the same type. Maybe this should
            #        instead raise an exception that's caught higher up...
            return None

        i0p, i1p = np.where(mask)[0][0], np.where(mask)[0][-1]
        f0 = (abs(hdu_waves[i0p] - fov_wmin + 0.5 * wdel) % wdel) / wdel    # blue edge
        f1 = (abs(hdu_waves[i1p] - fov_wmax - 0.5 * wdel) % wdel) / wdel    # red edge

        data = cubehdu.data[i0p:i1p+1,
                            xy0p[1]:xy1p[1],
                            xy0p[0]:xy1p[0]]
        data[0, :, :] *= f0.to(1).round(11).value  # decompose dimensionless
        if i1p > i0p:
            data[-1, :, :] *= f1.to(1).round(11).value

        # w0, w1 : the closest cube wavelengths outside the fov edge wavelengths
        # fov_waves : the fov edge wavelengths
        # f0, f1 : the scaling factors for the blue and red edge cube slices
        #
        # w0, w1 = hdu_waves[i0p], hdu_waves[i1p]

        new_hdr.update({
            "NAXIS": 3,
            "NAXIS3": data.shape[0],
            "CRVAL3": hdu_waves[i0p].to(hdr["CUNIT3"]).round(11).value,
            "CRPIX3": 1,  # 1 because FITS...
            "CDELT3": hdr["CDELT3"],
            "CUNIT3": hdr["CUNIT3"],
            "CTYPE3": hdr["CTYPE3"],
            "BUNIT":  hdr["BUNIT"],
        })

        return new_hdr, data

    def _calc_area_factor(self, field):
        bg_solid_angle = u.Unit(field.header["SOLIDANG"]).to(u.arcsec**-2)
        # arcsec**2 * arcsec**-2
        area_factor = self.pixel_area * bg_solid_angle
        return area_factor

    def _make_spectrum_cubefields(self):
        """
        Find Cube fields.

        * collapse cube along spatial dimensions --> spectrum vector
        * convert vector to PHOTLAM
        * interpolate at waveset
        * yield scaled flux to be added to canvas flux
        """
        for field in self._get_cube_fields():
            fluxes = field.data.sum(axis=(1, 2))
            fov_waveset_fluxes = np.interp(self.waveset, field.waveset, fluxes)

            # .lower() is needed because astropy doesn't recognize PHOTLAM
            field_unit = field.header.get("BUNIT", PHOTLAM).lower()
            flux_scale_factor = u.Unit(field_unit).to(PHOTLAM)

            yield fov_waveset_fluxes * flux_scale_factor

    def _make_spectrum_imagefields(self):
        """
        Find Image fields.

        * sum image over both dimensions
        * evaluate spectum at waveset
        * yield spectrum multiply by sum to be added to canvas flux
        """
        for field in self._get_image_fields():
            weight = np.sum(field.data)  # Shouldn't that be 1 by convention?
            yield field.spectrum(self.waveset).value * weight

    def _make_spectrum_tablefields(self):
        """
        Find Table fields.

        * evaluate all spectra at waveset
        * for each unique ref, sum the weights
        * yield each spectrum * sum of weights to be added to canvas flux
        """
        for field in self._get_table_fields():
            refs = np.array(field["ref"])
            weights = np.array(field["weight"])
            # TODO: could do grouping of table with both columns??
            for ref in set(refs):
                weight = np.sum(weights, where=refs == ref)
                yield field.spectra[int(ref)](self.waveset).value * weight

    def _make_spectrum_backfields(self):
        for field in self._get_background_fields():
            weight = self._calc_area_factor(field)
            yield field.spectrum(self.waveset).value * weight

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
        # Start with zero flux no ensure correct array shape even if none of
        # the sub-functions yield anything.
        canvas_flux = sum(chain(
            self._make_spectrum_cubefields(),
            self._make_spectrum_imagefields(),
            self._make_spectrum_tablefields(),
            self._make_spectrum_backfields(),
        ), start=np.zeros_like(self.waveset.value))

        spectrum = SourceSpectrum(Empirical1D, points=self.waveset,
                                  lookup_table=canvas_flux)
        return spectrum

    def _make_image_cubefields(self):
        """
        Find Cube fields.

        * collapse cube along wavelength axis
        * rescale and reorient image
        * yield cube image  to be added to canvas image
        """
        for field in self._get_cube_fields():
            # cube_fields come in with units of photlam/arcsec2,
            # need to convert to ph/s
            # We need to the voxel volume (spectral and solid angle) for that.
            # ..todo: implement branch for use_photlam is True
            spectral_bin_width = (field.header["CDELT3"] *
                                  u.Unit(field.header["CUNIT3"])
                                  ).to(u.Angstrom)
            # First collapse to image, then convert units
            image = np.sum(field.data, axis=0) * PHOTLAM/u.arcsec**2
            image = (image * self._pixarea(field.header) * self.area *
                     spectral_bin_width).to(u.ph/u.s)
            yield fits.ImageHDU(data=image, header=field.header)

    def _make_image_imagefields(self, fov_waveset, bin_widths, use_photlam):
        """
        Find Image fields.

        * sum spectra between wavelength edges
        * multiply image by summed spectrum
        * yield image  to be added to canvas image
        """
        for field in self._get_image_fields():
            image = deepcopy(field.data)

            # TODO: Improve this...
            spec = (field.spectrum(fov_waveset) if use_photlam
                    else (field.spectrum(fov_waveset) *
                          bin_widths * self.area).to(u.ph / u.s))

            flux = spec.value.sum()
            image *= flux  # ph / s
            yield fits.ImageHDU(data=image, header=field.header)

    def _make_image_tablefields(self, fov_waveset, bin_widths, use_photlam):
        """
        Find Table fields.

        * sum spectra between wavelength edges
        * yield summed flux at x,y position to be added to canvas image
        """
        for field in self._get_table_fields():
            # x, y are ALWAYS in arcsec - crval is in deg
            xpix, ypix = imp_utils.val2pix(self.header,
                                           field["x"] / 3600,
                                           field["y"] / 3600)

            fluxes = {
                ref: spec(fov_waveset).value.sum() if use_photlam
                else (spec(fov_waveset) *
                      bin_widths * self.area).to(u.ph / u.s).value.sum()
                for ref, spec in field.spectra.items()
            }

            if self.sub_pixel:
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

    def _make_image_backfields(self, fov_waveset, bin_widths, use_photlam):
        for field in self._get_background_fields():

            # TODO: Improve this...
            spec = (field.spectrum(fov_waveset) if use_photlam
                    else (field.spectrum(fov_waveset) *
                          bin_widths * self.area).to(u.ph / u.s))

            flux = spec.value.sum()
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
        # Make waveset and canvas image
        fov_waveset = self.waveset
        bin_widths = np.diff(fov_waveset)       # u.um
        bin_widths = 0.5 * (np.r_[0, bin_widths] + np.r_[bin_widths, 0])

        canvas_image_hdu = fits.ImageHDU(
            data=np.zeros((self.header["NAXIS2"], self.header["NAXIS1"])),
            header=self.header)

        for tmp_hdu in chain(self._make_image_cubefields(),
                             self._make_image_imagefields(
                                 fov_waveset, bin_widths, use_photlam)):
            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(
                tmp_hdu,
                canvas_image_hdu,
                conserve_flux=True,
                spline_order=self.spline_order)

        for flux, weight, x, y in self._make_image_tablefields(
                fov_waveset, bin_widths, use_photlam):
            if self.sub_pixel:
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

        canvas_image_hdu.data = sum(
            self._make_image_backfields(fov_waveset, bin_widths, use_photlam),
            start=canvas_image_hdu.data)

        canvas_image_hdu.header["BUNIT"] = "ph s-1"
        return canvas_image_hdu  # [ph s-1]

    def _make_cube_cubefields(self, fov_waveset):
        """
        Find Cube fields.

        * rescale and reorient cubes
        * interp1d smaller cubes with waveset
        * yield cubes to be added to cavas cube
        """
        for field in self._get_cube_fields():
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Assumption is that ImageHDUs have units of PHOTLAM arcsec-2

            # ..todo: Deal with this bounds_error in a more elegant way
            field_interp = interp1d(field.waveset.to(u.um).value,
                                    field.data, axis=0, kind="linear",
                                    bounds_error=False, fill_value=0)

            field_data = field_interp(fov_waveset.value)

            # Pixel scale conversion
            field_data *= self._pixarea(field.header).value / self.pixel_area
            field_hdu = fits.ImageHDU(data=field_data, header=field.header)
            yield field_hdu

    def _make_cube_imagefields(self, fov_waveset, spline_order):
        """
        Find Image fields.

        * rescale and reorient images
        * evaluate spectra at waveset
        * expand image by spectra to 3D form
        * yield image cubes to be added to cavas cube
        """
        for field in self._get_image_fields():
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Assumption is that ImageHDUs have units of PHOTLAM arcsec-2
            # ImageHDUs have photons/second/pixel.
            canvas_image_hdu = fits.ImageHDU(
                data=np.zeros((self.header["NAXIS2"], self.header["NAXIS1"])),
                header=self.header)
            # FIX: Do not scale source data - make a copy first.
            bunit = u.Unit(field.header.get("BUNIT", ""))
            field_data = deepcopy(field.data)
            if unit_includes_per_physical_type(bunit, "solid angle"):
                # Field is in (PHOTLAM) / arcsec**2, need to scale by pixarea
                field_data *= self._pixarea(field.header).value
            field_data /= self.pixel_area
            field_hdu = fits.ImageHDU(data=field_data, header=field.header)

            canvas_image_hdu = imp_utils.add_imagehdu_to_imagehdu(
                field_hdu,
                canvas_image_hdu,
                spline_order=spline_order)
            spec = field.spectrum(fov_waveset)

            # 2D * 1D -> 3D
            field_cube = canvas_image_hdu.data[None, :, :] * spec[:, None, None]
            yield field_cube.value

    def _make_cube_tablefields(self, fov_waveset):
        """
        Find Table fields.

        * evaluate spectra at waveset
        * yield spectrum at x,y position to be added cavas to cube
        """
        for field in self._get_table_fields():
            # Cube should be in PHOTLAM arcsec-2 for SpectralTrace mapping
            # Point sources are in PHOTLAM per pixel
            # Point sources need to be scaled up by inverse pixel_area
            for row in field:
                xsky, ysky = row["x"] / 3600, row["y"] / 3600
                # x, y are ALWAYS in arcsec - crval is in deg
                # TODO: Change this to some proper WCS function!
                xpix, ypix = imp_utils.val2pix(self.header, xsky, ysky)
                flux = field.spectra[row["ref"]](fov_waveset)
                flux_vector = flux.value * row["weight"] / self.pixel_area

                if self.sub_pixel:
                    xs, ys, fracs = imp_utils.sub_pixel_fractions(xpix, ypix)
                    for i, j, k in zip(xs, ys, fracs):
                        yield flux_vector * k, i, j
                else:
                    yield flux_vector, int(xpix), int(ypix)

    def _make_cube_backfields(self, fov_waveset):
        for field in self._get_background_fields():
            # FIXME: This assumes that SOLIDANG == arcsec-2, which is usually
            #        True, but doesn't have to be. Maybe solve via BUNIT?
            #        Remember, cube output needs PHOTLAM / arcsec**2 !
            #        So if SOLIDANG or BUNIT or whatever is not in arcsec-2,
            #        the spectrum shoule be scaled accordingly!
            spec = field.spectrum(fov_waveset)
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

        ``PHOTLAM = ph / (cm2 * s * AA)``.
        Original source fields are in units of:

        - tables: (PHOTLAM in spectrum)
        - images: arcsec-2 (PHOTLAM in spectrum)
        - cubes: PHOTLAM arcsec-2

        .. warning:: Input Images and Cubes should have units of PHOTLAM arcsec-2

        Returns
        -------
        canvas_cube_hdu : fits.ImageHDU
            [ph s-1 um-1 arcsec-2]      # as needed by SpectralTrace

        """
        # 1. Make waveset and canvas cube (area, bin_width are applied at end)
        # TODO: Why is this not self.waveset? What's different?
        # -> For non-cube input but cube output (e.g. 2D image in spec mode),
        #    waverange needs to be resampled (I guess) to spectral_bin_width,
        #    but self.waverange can only access the fields.
        # -> Perhaps change that??
        wave_unit = u.Unit(from_currsys("!SIM.spectral.wave_unit", self.cmds))
        fov_waveset = np.arange(
            self.meta["wave_min"].value, self.meta["wave_max"].value,
            from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)) * wave_unit

        # Note: There is an edge case where some floating point madness
        #       results in the arange ticking over by another bin by just
        #       .000000000005 (yeah...), which creates and off-by-one error
        #       further down the line.
        # TODO: Find a better way to solve this, perhaps with linspace...
        n_wave = int(
            (self.meta["wave_max"].value - self.meta["wave_min"].value) /
            from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)
        )
        fov_waveset = fov_waveset[:n_wave].to(u.um)

        # TODO: what's with this code??
        # fov_waveset = self.waveset
        # wave_bin_n = len(fov_waveset)
        # if "lin" in self.meta["wave_bin_type"]:
        #     fov_waveset = np.linspace(wave_min, wave_max, wave_bin_n)
        # elif "log" in self.meta["wave_bin_type"]:
        #     wmin, wmax = wave_min.to(u.um).value, wave_max.to(u.um).value
        #     fov_waveset = np.logspace(wmin, wmax, wave_bin_n)

        # make canvas cube based on waveset of largest cube and NAXIS1,2
        # from fov.header
        canvas_cube_hdu = fits.ImageHDU(
            data=np.zeros((len(fov_waveset),
                           self.header["NAXIS2"],
                           self.header["NAXIS1"])),
            header=self.header)
        # set BUNIT initially to PHOTLAM / arcsec**2
        canvas_cube_hdu.header["BUNIT"] = "ph cm-2 s-1 AA-1 arcsec-2"

        canvas_cube_hdu.header.update({
            "CDELT3": np.diff(fov_waveset[:2])[0].to_value(u.um),
            "CRVAL3": fov_waveset[0].value,
            "CRPIX3": 1,
            "CUNIT3": "um",
            "CTYPE3": "WAVE",
        })
        # TODO: Add the log wavelength keyword here, if log scale is needed

        if (self.detector_header is not None and
                self.detector_header["NAXIS"] == 3):  # cube ifu mode
            # TODO: Find out why not always do that, probably related to this
            #       "!INST.decouple_detector_from_sky_headers" thingy
            new_det_wcs = WCS(self.detector_header, key="D")
            canvas_cube_hdu.header.update(new_det_wcs.to_header())

            assert canvas_cube_hdu.header["NAXIS3"] == self.detector_header["NAXIS3"]

        for field_hdu in self._make_cube_cubefields(fov_waveset):
            canvas_cube_hdu = imp_utils.add_imagehdu_to_imagehdu(
                field_hdu,
                canvas_cube_hdu,
                spline_order=self.spline_order)

        canvas_cube_hdu.data = sum(self._make_cube_imagefields(
            fov_waveset, self.spline_order),
            start=canvas_cube_hdu.data)

        for flux, x, y in self._make_cube_tablefields(fov_waveset):
            # To prevent adding array values in this manner.
            assert not isinstance(x, Iterable), "x should be integer"
            canvas_cube_hdu.data[:, y, x] += flux

        canvas_cube_hdu.data = sum(self._make_cube_backfields(fov_waveset),
                                   start=canvas_cube_hdu.data)

        # TODO: what's with this code??
        # 6. Convert from PHOTLAM to ph/s/voxel
        #    PHOTLAM = ph/s/cm-2/AA
        #    area = m2, fov_waveset = um
        # SpectralTrace wants ph/s/um/arcsec2 --> get rid of m2, leave um
        canvas_cube_hdu.data *= self.area.to(u.cm ** 2).value
        canvas_cube_hdu.data *= 1e4       # ph/s/AA/arcsec2 --> ph/s/um/arcsec2
        canvas_cube_hdu.header["BUNIT"] = "ph s-1 um-1 arcsec-2"

        # TODO: what's with this code??
        # bin_widths = np.diff(fov_waveset).to(u.AA).value
        # bin_widths = 0.5 * (np.r_[0, bin_widths] + np.r_[bin_widths, 0])
        # canvas_cube_hdu.data *= bin_widths[:, None, None]

        return canvas_cube_hdu      # [ph s-1 um-1 (arcsec-2)]

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

    def get_corners(self, new_unit: str = None):
        """Return sky footprint, image plane footprint."""
        sky_corners = imp_utils.calc_footprint(self.header, new_unit=new_unit)
        imp_corners = imp_utils.calc_footprint(self.header, "D")
        return sky_corners, imp_corners

    @property
    def waverange(self):
        """Return wavelength range in um [wave_min, wave_max]."""
        if self._waverange is None:
            wave_min = self.meta["wave_min"]
            wave_max = self.meta["wave_max"]
            self._waverange = [wave_min, wave_max] << u.um
        return self._waverange

    @property
    def wavelength(self):
        """Return central wavelength in um."""
        if self._wavelength is None:
            self._wavelength = self.waverange.mean() << u.um
        return self._wavelength

    @property
    def waveset(self):
        """Return a wavelength vector in um."""
        if cube_fields := self._get_cube_fields():
            naxis3_max = np.argmax([cube.header["NAXIS3"]
                                    for cube in cube_fields])
            _waveset = cube_fields[naxis3_max].waveset
        elif specfields := self._get_fields(SpectrumSourceField):
            _waveset = np.concatenate([
                spec.waveset.to(u.um)
                for field in specfields
                for spec in field.spectra.values()
            ])
        else:
            _waveset = self.waverange << u.um

        # TODO: tie the round to a global precision setting (this is defined
        #       better in another TODO somewhere...)
        # The rounding is necessary to not end up with things like:
        #   0.7 um
        #   0.7000000000000001 um
        #   0.7000000000000002 um
        # and yes, that actually happend...
        _waveset = np.unique(_waveset.round(10))
        return _waveset

    @property
    def area(self) -> float:
        """Return meta["area"] from cmds."""
        # TODO: Make this an attribute once meta resolving is implemented.
        return from_currsys(self.meta["area"], self.cmds) << u.m**2

    @property
    def sub_pixel(self) -> bool:
        """Return meta["sub_pixel"] from cmds."""
        # TODO: Make this an attribute once meta resolving is implemented.
        return from_currsys(self.meta["sub_pixel"], self.cmds)

    @property
    def spline_order(self) -> int:
        """Return "!SIM.computing.spline_order" from cmds."""
        # TODO: Make this an attribute once meta resolving is implemented.
        return from_currsys("!SIM.computing.spline_order", self.cmds)

    def _get_fields(self, subclass) -> list[SourceField]:
        """Return list of fields of specific SourceField subclass."""
        return [field for field in self.fields if isinstance(field, subclass)]

    def _get_cube_fields(self) -> list[CubeSourceField]:
        """Return list of CubeSourceFields."""
        return self._get_fields(CubeSourceField)

    def _get_image_fields(self) -> list[ImageSourceField]:
        """Return list of ImageSourceFields."""
        return self._get_fields(ImageSourceField)

    def _get_table_fields(self) -> list[TableSourceField]:
        """Return list of TableSourceFields."""
        return self._get_fields(TableSourceField)

    def _get_background_fields(self) -> list[BackgroundSourceField]:
        """Return list of BackgroundSourceFields."""
        return self._get_fields(BackgroundSourceField)

    @property
    def spectra(self) -> dict[int, SourceSpectrum]:
        """Return a collection of all fields' spectra.

        .. deprecated:: 0.10.0

           Use individual fields' spectra instead.
        """
        warn("Deprecated since v0.10.0.",
             DeprecationWarning, stacklevel=2)
        specs = {
            ref: spec
            for field in [
                fld for fld in self.fields
                if isinstance(fld, SpectrumSourceField)
            ]
            for ref, spec in field.spectra.items()
        }
        return specs

    def _ensure_deg_header(self):
        cunit = u.Unit(self.header["CUNIT1"].lower())
        convf = cunit.to(u.deg)
        self.header["CDELT1"] *= convf
        self.header["CDELT2"] *= convf
        self.header["CRVAL1"] *= convf
        self.header["CRVAL2"] *= convf
        self.header["CUNIT1"] = "deg"
        self.header["CUNIT2"] = "deg"

    def plot(self, axes, units: str = "arcsec") -> None:
        """Plot FOV footprint."""
        outline = np.array(list(close_loop(self.get_corners(units)[0])))
        axes.plot(*outline.T, label=f"FOV id: {self.meta['id']}")

    def __repr__(self) -> str:
        """Return repr(self)."""
        msg = (f"{self.__class__.__name__}({self.header!r}, "
               f"{self.waverange!r}, {self.detector_header!r}, "
               f"**{self.meta!r})")
        return msg

    def __str__(self) -> str:
        """Return str(self)."""
        msg = (f"FOV id: {self.meta['id']}, with dimensions "
               f"({self.header['NAXIS1']}, {self.header['NAXIS2']})\n"
               f"Sky centre: ({self.header['CRVAL1']}, "
               f"{self.header['CRVAL2']})\n"
               f"Image centre: ({self.header['CRVAL1D']}, "
               f"{self.header['CRVAL2D']})\n"
               f"Wavelength range: {self.waverange!s}\n")
        return msg


def replace_nans(field, cmds) -> None:
    """Replace NaN values in image field."""
    try:
        fill_value = cmds["!SIM.computing.nan_fill_value"]
    except TypeError:
        # This occurs in testing, not sure why
        fill_value = 0.
    if np.ma.fix_invalid(field.data, copy=False,
                         fill_value=fill_value).mask.any():
        logger.warning("Source contained NaN values, which "
                       "were replaced by %f.", fill_value)
        logger.info(
            "The replacement value for NaNs in sources can be "
            "set in '!SIM.computing.nan_fill_value'.")


def extract_range_from_spectrum(spectrum, waverange):
    assert isinstance(spectrum, SourceSpectrum), (
        f"spectrum must be of type synphot.SourceSpectrum: {type(spectrum)}")

    wave_min, wave_max = quantify(waverange, u.um).to(u.AA).value
    spec_waveset = spectrum.waveset.to(u.AA).value
    mask = (spec_waveset > wave_min) * (spec_waveset < wave_max)

    # FIXME: Why did I comment this out in 2023? Seems useful to have...
    # if sum(mask) == 0:
    #     logger.info(
    #         "Waverange does not overlap with Spectrum waveset: %s <> %s for "
    #         "spectrum %s", [wave_min, wave_max], spec_waveset, spectrum)
    if wave_min < min(spec_waveset) or wave_max > max(spec_waveset):
        logger.info(("Waverange only partially overlaps with Spectrum waveset: "
                     "%s <> %s for spectrum %s"),
                     [wave_min, wave_max], spec_waveset, spectrum)

    wave = np.r_[wave_min, spec_waveset[mask], wave_max]
    flux = spectrum(wave)

    new_spectrum = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)
    new_spectrum.meta.update(spectrum.meta)

    return new_spectrum
