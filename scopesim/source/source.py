"""
# old functionality to implement:
# - provide x, y, lam, spectra, weight, ref
# - overridden + : number, Source, SourceSpectrum
# - overridden * : number, SpectralElement
# - write to and read from file
# - shift all fields
# - rotate around the centre
# - photons_in_range returns the photons per spectrum in a wavelength range
# - image_in_range returns an image of the source for a wavelength range
#
# old functionality which will be removed:
# - project_onto_chip
# - apply_optical_train
#
# old structure --> new structure:
# - all data held in 6 arrays
# --> new dicts for fields, spectrum
#       field can be a Table or an ImageHDU
#       spectrum is a SourceSpectrum
#
# Use cases:
# image + spectrum
# images + spectra
# table + spectrum
# table + spectra
#
# table columns = x, y, spec_id, weight
# table meta keywords = x_unit, y_unit
#
# image header keywords = WCS, SPEC_ID, WEIGHT
# [WCS = CRPIXn, CRVALn = (0,0), CTYPEn, CDn_m, NAXISn, CUNITn
"""

# import pickle
from copy import deepcopy
from pathlib import Path
from typing import Union, TextIO
from io import StringIO

import numpy as np

from astropy.table import Table, Column
from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

from synphot import SpectralElement, SourceSpectrum

from ..optics.image_plane import ImagePlane
from ..optics import image_plane_utils as imp_utils
from .source_utils import validate_source_input, convert_to_list_of_spectra, \
    photons_in_range
from . import source_templates as src_tmp

from ..base_classes import SourceBase
from ..utils import (find_file, is_fits, get_fits_type, quantify,
                     quantity_from_table, convert_table_comments_to_dict,
                     close_loop, figure_factory, get_logger)


logger = get_logger(__name__)


class Source(SourceBase):
    """
    Create a source object from a file or from arrays.

    A Source object must consist of a spatial and a spectral description
    of the on-sky source. Many sources can be added together and kept
    in memory as a single Source object.

    The spatial descriptions are kept in the ``<Source>.fields`` list,
    while the spectral descriptions are in the ``<Source>.spectra`` list.

    The spatial description can be built from any combination of:

    * a list of arrays (like in SimCADO >v0.5)
    * astropy Table objects
    * astropy ImageHDU objects
    * on disk FITS files
    * on disk ASCII tables

    The spectral descriptions can be passed as either ``synphot.SourceSpectrum``
    objects, or a set of two equal length arrays for wavelength and flux.

    .. hint:: Initialisation parameter combinations include:

        New ScopeSim-style input
        - ``table=<astropy.Table>, spectra=<list of synphot.SourceSpectrum>``
        - ``table=<astropy.Table>, lam=<array>, spectra=<list of array>``
        - ``image_hdu=<fits.ImageHDU>, spectra=<list of synphot.SourceSpectrum>``
        - ``image_hdu=<fits.ImageHDU>, lam=<array>, spectra=<list of array>``
        - ``image_hdu=<fits.ImageHDU>, flux=<astropy.Quantity>``

        Old SimCADO-style input
        - ``x=<array>, y=<array>, ref=<array>, spectra=<list of synphot.SourceSpectrum>``
        - ``x=<array>, y=<array>, ref=<array>, spectra=<list of array>, lam=<array>``
        - ``x=<array>, y=<array>, ref=<array>, weight=<array>, spectra=<list of array>, lam=<array>``

        More details on the content of these combinations can be found in the
        use-case documentation.

    Parameters
    ----------
    filename : str

    lam : np.array
        [um] Wavelength bins of length (m)
    spectra : list of synphot.SourceSpectra
        [ph/s/cm2/AA]
    x, y : np.array
        [arcsec] coordinates of where the emitting files are relative to the
        centre of the field of view
    ref : np.array
        the index for .spectra which connects a position (x, y) to a spectrum
        ``flux(x[i], y[i]) = spectra[ref[i]] * weight[i]``
    weight : np.array
        A weighting to scale the relevant spectrum for each position
    table : astropy.Table
    image_hdu : fits.ImageHDU
        [arcsec-2] The .data array is simply a map of weights for the assiciated
        spectrum referenced by .header["SPEC_REF].
        Surface brightness values are assumed to be per arcsec2
    flux : astropy.Quantity
        [u.mag, u.ABmag, u.Jy] Flux values are converted to a reference spectrum
        that is referenced by image_hdu.header["SPEC_REF"].
        flux can only be used in conjuction with image_hdu

    Attributes
    ----------
    fields : list
        The spatial distribution of the on-sky source, either as
        ``fits.ImageHDU`` or ``astropy.Table`` objects
    spectra : list of ``synphot.SourceSpectrum`` objects
        List of spectra associated with the fields
    meta : dict
        Dictionary of extra information about the source

    See Also
    --------
    synphot : ``https://synphot.readthedocs.io/en/latest/``

    """

    def __init__(self, filename=None, cube=None, ext=0,
                 lam=None, spectra=None, x=None, y=None, ref=None, weight=None,
                 table=None, image_hdu=None, flux=None, **kwargs):

        self.meta = {}
        self.meta.update(kwargs)
        # ._meta_dicts contains a meta for each of the .fields. It is primarily
        # used to set proper FITS header keywords for each field so the source
        # can be reconstructed from the FITS headers.
        self._meta_dicts = [self.meta]

        self.fields: list[Union[Table, fits.ImageHDU]] = []
        self.spectra: dict[int, SourceSpectrum] = {}

        self.bandpass = None

        validate_source_input(lam=lam, x=x, y=y, ref=ref, weight=weight,
                              spectra=spectra, table=table, cube=cube,
                              ext=ext, image_hdu=image_hdu, flux=flux,
                              filename=filename)

        if spectra is not None:
            spectra = {i: spec for i, spec in
                       enumerate(convert_to_list_of_spectra(spectra, lam))}

        if filename is not None and spectra is not None:
            self._from_file(filename, spectra, flux)

        elif cube is not None:
            self._from_cube(cube=cube, ext=ext)

        elif table is not None and spectra is not None:
            self._from_table(table, spectra)

        elif image_hdu is not None and spectra is not None:
            self._from_imagehdu_and_spectra(image_hdu, spectra)

        elif image_hdu is not None and flux is not None:
            self._from_imagehdu_and_flux(image_hdu, flux)

        elif image_hdu is not None and flux is None and spectra is None:
            if image_hdu.header.get("BUNIT") is not None:
                self._from_imagehdu_only(image_hdu)
            else:
                msg = ("image_hdu must be accompanied by either spectra or "
                       f"flux:\n spectra: {spectra}, flux: {flux}")
                raise ValueError(msg)

        elif x is not None and y is not None and \
                ref is not None and spectra is not None:
            self._from_arrays(x, y, ref, weight, spectra)

    def _from_file(self, filename, spectra, flux):
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        filename = find_file(filename)

        if not is_fits(filename):
            tbl = ioascii.read(filename)
            tbl.meta.update(convert_table_comments_to_dict(tbl))
            self._from_table(tbl, spectra)
            return

        fits_type = get_fits_type(filename)
        data = fits.getdata(filename)
        hdr = fits.getheader(filename)
        hdr["FILENAME"] = Path(filename).name
        if fits_type == "image":
            image = fits.ImageHDU(data=data, header=hdr)
            if spectra is not None:
                self._from_imagehdu_and_spectra(image, spectra)
            elif flux is not None:
                self._from_imagehdu_and_flux(image, flux)
            else:
                self._from_imagehdu_only(image)
        elif fits_type == "bintable":
            hdr1 = fits.getheader(filename, 1)
            hdr.update(hdr1)
            tbl = Table(data, meta=dict(hdr))
            tbl.meta.update(convert_table_comments_to_dict(tbl))
            self._from_table(tbl, spectra)

    def _from_table(self, tbl: Table, spectra: dict[int, SourceSpectrum]):
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        if "weight" not in tbl.colnames:
            tbl.add_column(Column(name="weight", data=np.ones(len(tbl))))
        # tbl["ref"] += len(self.spectra)
        # self.fields.append(tbl)
        # self.spectra += spectra

        self.fields = [tbl]
        self.spectra = spectra

    def _from_imagehdu_and_spectra(self, image_hdu, spectra):
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        if not image_hdu.header.get("BG_SRC"):
            pass
            # FIXME: This caused more problems than it solved!
            #        Find out if there's a good reason to mess with this,
            #        otherwise just remove...

            # image_hdu.header["CRVAL1"] = 0
            # image_hdu.header["CRVAL2"] = 0
            # image_hdu.header["CRPIX1"] = image_hdu.header["NAXIS1"] / 2
            # image_hdu.header["CRPIX2"] = image_hdu.header["NAXIS2"] / 2
            # #image_hdu.header["CRPIX1"] = (image_hdu.header["NAXIS1"] + 1) / 2
            # #image_hdu.header["CRPIX2"] = (image_hdu.header["NAXIS2"] + 1) / 2
            # # .. todo:: find where the actual problem is with negative CDELTs
            # # .. todo:: --> abs(pixel_scale) in header_from_list_of_xy
            # if image_hdu.header["CDELT1"] < 0:
            #     image_hdu.header["CDELT1"] *= -1
            #     image_hdu.data = image_hdu.data[:, ::-1]
            # if image_hdu.header["CDELT2"] < 0:
            #     image_hdu.header["CDELT2"] *= -1
            #     image_hdu.data = image_hdu.data[::-1, :]

        if isinstance(image_hdu, fits.PrimaryHDU):
            image_hdu = fits.ImageHDU(data=image_hdu.data,
                                      header=image_hdu.header)

        if not spectra:
            raise ValueError("No spectrum was provided.")

        # image_hdu.header["SPEC_REF"] = len(self.spectra)
        assert len(spectra) == 1, f"_from_imagehdu_and_spectra needs single spectrum, ref was {image_hdu.header.get('SPEC_REF')}"
        image_hdu.header["SPEC_REF"] = 0
        # self.spectra += spectra
        self.spectra = spectra

        for i in [1, 2]:
            # Do not test for CUNIT or CDELT so that it throws an exception
            unit = u.Unit(image_hdu.header[f"CUNIT{i}"].lower())
            val = float(image_hdu.header[f"CDELT{i}"])
            image_hdu.header[f"CUNIT{i}"] = "deg"
            image_hdu.header[f"CDELT{i}"] = val * unit.to(u.deg)

        # self.fields.append(image_hdu)
        self.fields = [image_hdu]

    def _from_imagehdu_and_flux(self, image_hdu, flux):
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        if isinstance(flux, u.Unit):
            flux = 1 * flux

        spec_template = src_tmp.vega_spectrum
        if isinstance(flux, u.Quantity):
            if flux.unit.physical_type == "spectral flux density":  # ABmag and Jy
                spec_template = src_tmp.ab_spectrum
                flux = flux.to(u.ABmag)
            flux = flux.value
        spectra = {0: spec_template(flux)}
        self._from_imagehdu_and_spectra(image_hdu, spectra)

    def _from_imagehdu_only(self, image_hdu):
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        bunit = image_hdu.header.get("BUNIT")
        try:
            bunit = u.Unit(bunit)
        except ValueError:
            logger.error(
                "Astropy cannot parse BUNIT [%s].\n"
                "You can bypass this check by passing an astropy Unit to the "
                "flux parameter:\n"
                ">>> Source(image_hdu=..., flux=u.Unit(bunit), ...)", bunit)

        value = 0 if bunit in [u.mag, u.ABmag] else 1
        self._from_imagehdu_and_flux(image_hdu, value * bunit)

    def _from_arrays(self, x, y, ref, weight,
                     spectra: dict[int, SourceSpectrum]):
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        if weight is None:
            weight = np.ones(len(x))

        x = quantify(x, u.arcsec)
        y = quantify(y, u.arcsec)
        # tbl = Table(names=["x", "y", "ref", "weight"],
        #             data=[x, y, np.array(ref) + len(self.spectra), weight])
        tbl = Table(names=["x", "y", "ref", "weight"],
                    data=[x, y, ref, weight])
        tbl.meta["x_unit"] = "arcsec"
        tbl.meta["y_unit"] = "arcsec"

        # self.fields.append(tbl)
        # self.spectra += spectra
        self.fields = [tbl]
        self.spectra = spectra

    def _from_cube(self, cube, ext=0):
        """

        Parameters
        ----------
        cube: a file, HDUList or a PrimaryHDU object containing the cube
        ext: int
            the extension where the cube is located if applicable.

        """
        if self.spectra or self.fields:
            raise ValueError("Constructor method must act on empty instance!")

        if isinstance(cube, fits.HDUList):
            data = cube[ext].data
            header = cube[ext].header
            wcs = WCS(cube[ext], fobj=cube)
        elif isinstance(cube, (fits.PrimaryHDU, fits.ImageHDU)):
            data = cube.data
            header = cube.header
            wcs = WCS(cube)
        else:
            with fits.open(cube) as hdul:
                data = hdul[ext].data
                header = hdul[ext].header
                header["FILENAME"] = Path(cube).name
                wcs = WCS(cube)

        try:
            bunit = header["BUNIT"]
            u.Unit(bunit)
        except KeyError:
            bunit = "erg / (s cm2 arcsec2)"
            logger.warning(
                "Keyword \"BUNIT\" not found, setting to %s by default", bunit)
        except ValueError as error:
            logger.error("\"BUNIT\" keyword is malformed: %s", error)
            raise

        # Compute the wavelength vector. This will be attached to the cube_hdu
        # as a new `wave` attribute.  This is not optimal coding practice.
        wave = wcs.all_pix2world(header["CRPIX1"], header["CRPIX2"],
                                 np.arange(data.shape[0]), 0)[-1]

        wave = (wave * u.Unit(wcs.wcs.cunit[-1])).to(u.um,
                                                     equivalencies=u.spectral())

        # WCS keywords must be updated because astropy.wcs converts wavelengths to 'm'
        header.update(wcs.to_header())

        target_cube = data
        target_hdr = header.copy()
        target_hdr["BUNIT"] = bunit

        cube_hdu = fits.ImageHDU(data=target_cube, header=target_hdr)
        cube_hdu.wave = wave          # ..todo: review wave attribute, bad practice

        # self.fields.append(cube_hdu)
        self.fields = [cube_hdu]

    @property
    def table_fields(self):
        """List of fields that are defined through tables"""
        fields = [field for field in self.fields if isinstance(field, Table)]
        return fields

    @property
    def image_fields(self):
        """List of fields that are defined through two-dimensional images"""
        fields = [field for field in self.fields if
                  isinstance(field, fits.ImageHDU) and
                  field.header["NAXIS"] == 2]
        return fields

    @property
    def cube_fields(self):
        """List of fields that are defined through three-dimensional cubes"""
        fields = [field for field in self.fields if
                  isinstance(field, fits.ImageHDU) and
                  field.header["NAXIS"] == 3]
        return fields

    # ..todo: rewrite this method
    def image_in_range(self, wave_min, wave_max, pixel_scale=1*u.arcsec,
                       layers=None, area=None, spline_order=1, sub_pixel=False):
        if layers is None:
            layers = range(len(self.fields))
        fields = [self.fields[ii] for ii in layers]

        hdr = imp_utils.get_canvas_header(fields, pixel_scale=pixel_scale)
        im_plane = ImagePlane(hdr)

        for field in fields:
            if isinstance(field, Table):
                fluxes = self.photons_in_range(wave_min, wave_max, area,
                                               field["ref"]) * field["weight"]
                x = quantity_from_table("x", field, u.arcsec)
                y = quantity_from_table("y", field, u.arcsec)
                tbl = Table(names=["x", "y", "flux"], data=[x, y, fluxes])
                tbl.meta.update(field.meta)
                hdu_or_table = tbl

            elif isinstance(field, fits.ImageHDU):
                if field.header["SPEC_REF"] != "":
                    ref = [field.header["SPEC_REF"]]
                    flux = self.photons_in_range(wave_min, wave_max, area, ref)
                    # [ph s-1] or [ph s-1 m-2] come out of photons_in_range

                # ## ..todo: CATCH UNITS HERE. DEAL WITH THEM PROPERLY
                # Currently assuming that all images are scaled appropriately
                # and that they have SPEC_REF

                # else:
                #     field = scale_imagehdu(field, area=area,
                #                            solid_angle=pixel_scale**2,
                #                            waverange=(wave_min, wave_max))
                #     # [ph s-1] or [ph s-1 m-2] come out of photons_in_range
                #     flux = 1

                image = field.data * flux
                hdu = fits.ImageHDU(header=field.header, data=image)
                hdu_or_table = hdu
            else:
                continue

            im_plane.add(hdu_or_table, sub_pixel=sub_pixel,
                         spline_order=spline_order)

        return im_plane

    def photons_in_range(self, wave_min, wave_max, area=None, indexes=None):
        """

        Parameters
        ----------
        wave_min : float, u.Quantity
            [um]
        wave_max : float, u.Quantity
            [um]
        area : float, u.Quantity, optional
            [m2]
        indexes : list of integers, optional

        Returns
        -------
        counts : u.Quantity list
            [ph / s / m2] if area is None
            [ph / s] if area is passed

        """
        if indexes is None:
            # indexes = range(len(self.spectra))
            indexes = self.spectra.keys()

        spectra = [self.spectra[ii] for ii in indexes]
        counts = photons_in_range(spectra, wave_min, wave_max, area=area,
                                  bandpass=self.bandpass)
        return counts

    def fluxes(self, wave_min, wave_max, **kwargs):
        return self.photons_in_range(wave_min, wave_max, **kwargs)

    def image(self, wave_min, wave_max, **kwargs):
        return self.image_in_range(wave_min, wave_max, **kwargs)

    # @classmethod
    # def load(cls, filename):
    #     """Load :class:'.Source' object from filename"""
    #     with open(filename, "rb") as fp1:
    #         src = pickle.load(fp1)
    #     return src

    # def dump(self, filename):
    #     """Save to filename as a pickle"""
    #     with open(filename, "wb") as fp1:
    #         pickle.dump(self, fp1)

    def shift(self, dx: float = 0, dy: float = 0, layers=None) -> None:
        """
        Shift the position of one or more fields w.r.t. the optical axis.

        Parameters
        ----------
        dx, dy : float
            [arcsec]
        layers : list of ints
            which .fields entries to shift

        """
        if layers is None:
            layers = np.arange(len(self.fields))

        for ii in layers:
            if isinstance(self.fields[ii], Table):
                x = quantity_from_table("x", self.fields[ii], u.arcsec)
                x += quantify(dx, u.arcsec)
                self.fields[ii]["x"] = x

                y = quantity_from_table("y", self.fields[ii], u.arcsec)
                y += quantify(dy, u.arcsec)
                self.fields[ii]["y"] = y
            elif isinstance(self.fields[ii], (fits.ImageHDU, fits.PrimaryHDU)):
                dx = quantify(dx, u.arcsec).to(u.deg)
                dy = quantify(dy, u.arcsec).to(u.deg)
                self.fields[ii].header["CRVAL1"] += dx.value
                self.fields[ii].header["CRVAL2"] += dy.value

    def rotate(self, angle, offset=None, layers=None):
        pass

    def add_bandpass(self, bandpass):
        if not isinstance(bandpass, SpectralElement):
            raise ValueError("type(bandpass) must be synphot.SpectralElement")

        self.bandpass = bandpass

    def plot(self):
        """
        Plot the location of source components.

        Source components instantiated from 2d or 3d ImageHDUs are represented
        by their spatial footprint. Source components instantiated from tables
        are shown as points.
        """
        _, axes = figure_factory()

        colours = "rgbcymk" * (len(self.fields) // 7 + 1)
        for col, field in zip(colours, self.fields):
            if isinstance(field, Table):
                axes.plot(field["x"], field["y"], col+".")
            elif isinstance(field, (fits.ImageHDU, fits.PrimaryHDU)):
                xypts = imp_utils.calc_footprint(field.header)
                convf = u.Unit(field.header["CUNIT1"]).to(u.arcsec)
                outline = np.array(list(close_loop(xypts))) * convf
                axes.plot(outline[:, 0], outline[:, 1], col)
                axes.set_xlabel("x [arcsec]")
                axes.set_ylabel("y [arcsec]")
        axes.set_aspect("equal")
        return axes

    def make_copy(self):
        new_source = Source()
        new_source.meta = deepcopy(self.meta)
        new_source._meta_dicts = deepcopy(self._meta_dicts)
        new_source.spectra = deepcopy(self.spectra)
        for field in self.fields:
            if isinstance(field, (fits.ImageHDU, fits.PrimaryHDU)) \
                    and field._file is not None:  # and field._data_loaded is False:
                new_source.fields.append(field)
            else:
                new_source.fields.append(deepcopy(field))

        return new_source

    def append(self, source_to_add):
        if not isinstance(source_to_add, Source):
            raise ValueError(f"Cannot add {type(source_to_add)} object to Source object")

        new_source = source_to_add.make_copy()
        # If there is no field yet, then self._meta_dicts contains a
        # reference to self.meta, which is empty. This ensures that both are
        # updated at the same time. However, it is important that the fields
        # and _meta_dicts match when appending sources.
        if len(self.fields) == 0:
            assert self._meta_dicts == [{}]
            self._meta_dicts = []

        specrefoffset = max(self.spectra.keys()) + 1 if self.spectra else 0
        for field in new_source.fields:
            if isinstance(field, Table):
                field["ref"] += specrefoffset
                self.fields.append(field)

            elif isinstance(field, (fits.ImageHDU, fits.PrimaryHDU)):
                if ("SPEC_REF" in field.header and
                        isinstance(field.header["SPEC_REF"], int)):
                    field.header["SPEC_REF"] += specrefoffset
                self.fields.append(field)

            self.spectra.update({k + specrefoffset: v for k, v in
                                 new_source.spectra.items()})
            self._meta_dicts += source_to_add._meta_dicts

    def __add__(self, new_source):
        self_copy = self.make_copy()
        self_copy.append(new_source)
        return self_copy

    def __radd__(self, new_source):
        return self.__add__(new_source)

    @staticmethod
    def _write_tablefield(fld: Table, stream: TextIO) -> None:
        stream.write(f"Table with {len(fld)} rows, referencing "
                     f"spectra {set(fld['ref'])}")

    @staticmethod
    def _write_imagefield(fld: Union[fits.ImageHDU, fits.PrimaryHDU],
                          stream: TextIO) -> None:
        im_size = fld.data.shape if fld.data is not None else "<empty>"
        stream.write(f"ImageHDU with size {im_size}, referencing "
                     f"spectrum {fld.header.get('SPEC_REF', '-')}")

    def __repr__(self) -> str:
        with StringIO() as str_stream:
            for ifld, fld in enumerate(self.fields):
                str_stream.write(f"[{ifld}]: ")
                if isinstance(fld, Table):
                    self._write_tablefield(fld, str_stream)
                elif isinstance(fld, (fits.ImageHDU, fits.PrimaryHDU)):
                    self._write_imagefield(fld, str_stream)
                str_stream.write("\n")
            return str_stream.getvalue()
