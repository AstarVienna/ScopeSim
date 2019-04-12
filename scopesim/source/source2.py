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

import pickle
import warnings
from copy import deepcopy
import numpy as np

from astropy.table import Table, Column
from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy import units as u

from synphot import SpectralElement

from ..optics import image_plane as opt_imp
from ..optics import image_plane_utils as imp_utils
from .source2_utils import validate_source_input, convert_to_list_of_spectra, \
    photons_in_range, scale_imagehdu

from .. import utils


class Source:
    """
    Create a source object from a file or from arrays
    
    A Source object must consist of a spatial and a spectral description
    of the on-sky source. Many sources can be added together and kept
    in memory as a single Source object.
    
    The spatial descriptions are kept in the ``<Source>.fields`` list,
    while the spectral descriptions are in ``<Source>.spectra`` list.
    
    The spatial description can be built from any combination of:
    
    * a list of arrays (like in SimCADO v0.5)
    * astropy Table objects
    * astropy ImageHDU objects
    * on disk FITS files
    * on disk ASCII tables
    
    while the spectral descriptions can be passed as either ``synphot.SourceSpectrum``
    objects, or a set of two equal length arrays for wavelength and flux.
   
    Parameters
    ----------
    spectra : array, or 
    filename : str
        
    or
    
    lam : np.array
        [um] Wavelength bins of length (m)
    spectra : np.array
        [ph/s/m2/bin] A (n, m) array with n spectra, each with m spectral bins
    x, y : np.array
        [arcsec] coordinates of where the emitting files are relative to the
        centre of the field of view
    ref : np.array
        the index for .spectra which connects a position (x, y) to a spectrum
        f(x[i], y[i]) = spectra[ref[i]] * weight[i]
    weight : np.array
        A weighting to scale the relevant spectrum for each position
    
    See Also
    --------
    ``synphot`` - https://synphot.readthedocs.io/en/latest/
    
    
    """

    def __init__(self, filename=None,
                 lam=None, spectra=None, x=None, y=None, ref=None, weight=None,
                 table=None, image_hdu=None, **kwargs):

        self.meta = {}
        self.meta.update(kwargs)
        
        self.fields = []
        self.spectra = []

        self.bandpass = None

        valid = validate_source_input(lam=lam, x=x, y=y, ref=ref, weight=weight,
                                      spectra=spectra, table=table,
                                      image_hdu=image_hdu, filename=filename)

        spectra = convert_to_list_of_spectra(spectra, lam)

        if filename is not None and spectra is not None:
            self._from_file(filename, spectra)

        elif table is not None and spectra is not None:
            self._from_table(table, spectra)

        elif image_hdu is not None and spectra is not None:
            self._from_imagehdu(image_hdu, spectra)

        elif x is not None and y is not None and \
                ref is not None and spectra is not None:
            self._from_arrays(x, y, ref, weight, spectra)

    def _from_file(self, filename, spectra):
        filename = utils.find_file(filename)

        if utils.is_fits(filename):
            fits_type = utils.get_fits_type(filename)
            data = fits.getdata(filename)
            hdr = fits.getheader(filename)
            if fits_type == "image":
                image = fits.ImageHDU(data=data, header=hdr)
                self._from_imagehdu(image, spectra)
            elif fits_type == "bintable":
                hdr1 = fits.getheader(filename, 1)
                hdr.update(hdr1)
                tbl = Table(data, meta=dict(hdr))
                tbl.meta.update(utils.convert_table_comments_to_dict(tbl))
                self._from_table(tbl, spectra)
        else:
            tbl = ioascii.read(filename)
            tbl.meta.update(utils.convert_table_comments_to_dict(tbl))
            self._from_table(tbl, spectra)

    def _from_table(self, tbl, spectra):
        if "weight" not in tbl.colnames:
            tbl.add_column(Column(name="weight", data=np.ones(len(tbl))))
        tbl["ref"] += len(self.spectra)
        self.fields += [tbl]
        self.spectra += spectra

    def _from_imagehdu(self, image_hdu, spectra):
        if spectra is not None and len(spectra) > 0:
            image_hdu.header["SPEC_REF"] = len(self.spectra)
            self.spectra += spectra
        else:
            image_hdu.header["SPEC_REF"] = ""
            warnings.warn("No spectrum was provided. SPEC_REF set to ''. "
                          "This could cause problems later")
            raise NotImplementedError

        self.fields += [image_hdu]

    def _from_arrays(self, x, y, ref, weight, spectra):
        if weight is None:
            weight = np.ones(len(x))

        x = utils.quantify(x, u.arcsec)
        y = utils.quantify(y, u.arcsec)
        tbl = Table(names=["x", "y", "ref", "weight"],
                    data=[x, y, np.array(ref) + len(self.spectra), weight])
        tbl.meta["x_unit"] = "arcsec"
        tbl.meta["y_unit"] = "arcsec"

        self.fields += [tbl]
        self.spectra += spectra

    def image_in_range(self, wave_min, wave_max, pixel_scale=1*u.arcsec,
                       layers=None, area=None, order=1, sub_pixel=False):
        if layers is None:
            layers = range(len(self.fields))
        fields = [self.fields[ii] for ii in layers]

        hdr = imp_utils.get_canvas_header(fields, pixel_scale=pixel_scale)
        im_plane = opt_imp.ImagePlane(hdr)

        for field in fields:
            if isinstance(field, Table):
                fluxes = self.photons_in_range(wave_min, wave_max, area,
                                               field["ref"]) * field["weight"]
                x = utils.quantity_from_table("x", field, u.arcsec)
                y = utils.quantity_from_table("y", field, u.arcsec)
                tbl = Table(names=["x", "y", "flux"], data=[x, y, fluxes])
                tbl.meta.update(field.meta)
                hdu_or_table = tbl

            elif isinstance(field, fits.ImageHDU):
                if field.header["SPEC_REF"] is not "":
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

            im_plane.add(hdu_or_table, sub_pixel=sub_pixel, order=order)

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
            indexes = range(len(self.spectra))

        spectra = [self.spectra[ii] for ii in indexes]
        counts = photons_in_range(spectra, wave_min, wave_max, area=area,
                                  bandpass=self.bandpass)
        return counts

    def fluxes(self, wave_min, wave_max, **kwargs):
        return self.photons_in_range(wave_min, wave_max, **kwargs)

    def image(self, wave_min, wave_max, **kwargs):
        return self.image_in_range(wave_min, wave_max, **kwargs)

    @classmethod
    def load(cls, filename):
        """Load :class:'.Source' object from filename"""
        with open(filename, 'rb') as fp1:
            src = pickle.load(fp1)
        return src

    def dump(self, filename):
        """Save to filename as a pickle"""
        with open(filename, 'wb') as fp1:
            pickle.dump(self, fp1)

    def shift(self, dx=0, dy=0, layers=None):

        if layers is None:
            layers = np.arange(len(self.fields))

        for ii in layers:
            if isinstance(self.fields[ii], Table):
                x = utils.quantity_from_table("x", u.arcsec)
                x += utils.quantify(dx, u.arcsec)
                self.fields[ii]["x"] = x

                y = utils.quantity_from_table("y", u.arcsec)
                y += utils.quantify(dy, u.arcsec)
                self.fields[ii]["y"] = y
            elif isinstance(self.fields[ii], fits.ImageHDU):
                dx = utils.quantify(dx, "arcsec").to(u.deg)
                dy = utils.quantify(dy, "arcsec").to(u.deg)
                self.fields[ii].header["CRVAL1"] += dx.value
                self.fields[ii].header["CRVAL2"] += dy.value

    def rotate(self, angle, offset=None, layers=None):
        pass

    def add_bandpass(self, bandpass):
        if not isinstance(bandpass, SpectralElement):
            raise ValueError("type(bandpass) must be synphot.SpectralElement")

        self.bandpass = bandpass

    def append(self, source_to_add):
        new_source = deepcopy(source_to_add)
        if isinstance(new_source, Source):
            for field in new_source.fields:
                if isinstance(field, Table):
                    field["ref"] += len(self.spectra)
                    self.fields += [field]

                elif isinstance(field, fits.ImageHDU):
                    if isinstance(field.header["SPEC_REF"], int):
                        field.header["SPEC_REF"] += len(self.spectra)
                    self.fields += [field]
                self.spectra += new_source.spectra
        else:
            raise ValueError("Cannot add {} object to Source object"
                             "".format(type(new_source)))

    def __add__(self, new_source):
        self_copy = deepcopy(self)
        self_copy.append(new_source)
        return self_copy

    def __radd__(self, new_source):
        return self.__add__(new_source)

    def __repr__(self):
        msg = ""
        for ii in range(len(self.fields)):
            if isinstance(self.fields[ii], Table):
                tbl_len = len(self.fields[ii])
                num_spec = set(self.fields[ii]["ref"])
                msg += "[{}]: Table with {} rows, referencing spectra {} \n" \
                       "".format(ii, tbl_len, num_spec)
            elif isinstance(self.fields[ii], fits.ImageHDU):
                im_size = self.fields[ii].data.shape
                num_spec = "-"
                if self.fields[ii].header["SPEC_REF"] is not "":
                    num_spec = self.fields[ii].header["SPEC_REF"]
                msg += "[{}]: ImageHDU with size {}, referencing spectrum {}" \
                       "\n".format(ii, im_size, num_spec)

        return msg
