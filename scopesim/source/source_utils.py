import warnings

import numpy as np
from astropy import wcs, units as u
from astropy.io import fits
from astropy.table import Table
import synphot
from synphot import units, SourceSpectrum, SpectralElement, Observation, Empirical1D
from synphot.models import Box1D
from synphot.units import convert_flux

from .. import utils


def validate_source_input(**kwargs):
    if "filename" in kwargs and kwargs["filename"] is not None:
        filename = kwargs["filename"]
        if utils.find_file(filename) is None:
            warnings.warn("filename was not found: {}".format(filename))

    if "image" in kwargs and kwargs["image_hdu"] is not None:
        image_hdu = kwargs["image_hdu"]
        if not isinstance(image_hdu, (fits.PrimaryHDU, fits.ImageHDU)):
            raise ValueError("image_hdu must be fits.HDU object with a WCS."
                             "type(image) == {}".format(type(image_hdu)))

        if len(wcs.find_all_wcs(image_hdu.header)) == 0:
            warnings.warn("image_hdu does not contain valid WCS. {}"
                          "".format(wcs.WCS(image_hdu)))

    if "table" in kwargs and kwargs["table"] is not None:
        tbl = kwargs["table"]
        if not isinstance(tbl, Table):
            raise ValueError("table must be an astropy.Table object:"
                             "{}".format(type(tbl)))

        if not np.all([col in tbl.colnames for col in ["x", "y", "ref"]]):
            raise ValueError("table must contain at least column names: "
                             "'x, y, ref': {}".format(tbl.colnames))

    return True


def convert_to_list_of_spectra(spectra, lam):
    spectra_list = []
    if isinstance(spectra, SourceSpectrum):
        spectra_list += [spectra]

    elif lam is None and\
            isinstance(spectra, (tuple, list)) and \
            isinstance(spectra[0], SourceSpectrum):
        spectra_list += spectra

    elif lam is not None and len(spectra.shape) == 1 and \
            isinstance(spectra, np.ndarray) and \
            isinstance(lam, np.ndarray):
        spec = SourceSpectrum(Empirical1D, points=lam, lookup_table=spectra)
        spectra_list += [spec]

    elif ((isinstance(spectra, np.ndarray) and
           len(spectra.shape) == 2) or
          (isinstance(spectra, (list, tuple)) and
           isinstance(spectra[0], np.ndarray))) and \
            isinstance(lam, np.ndarray):

        for sp in spectra:
            spec = SourceSpectrum(Empirical1D, points=lam, lookup_table=sp)
            spectra_list += [spec]

    return spectra_list


def get_vega_spectrum():
    """
    Copied from SimCADO
    Retrieve the Vega spectrum from stsci and return it in synphot format


    Notes
    -----
    To access wavelength and fluxes use::

        wave, flux = vega_sp._get_arrays(wavelengths=None)

    """
    location = "http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_008.fits"
    remote = synphot.specio.read_remote_spec(location, cache=True)
    header = remote[0]
    wave = remote[1]
    flux = remote[2]
    url = 'Vega from ' + location
    meta = {'header': header, 'expr': url}
    vega_sp = synphot.SourceSpectrum(Empirical1D, points=wave,
                                     lookup_table=flux, meta=meta)
    return vega_sp


def photons_in_range(spectra, wave_min, wave_max, area=None, bandpass=None):
    """

    Parameters
    ----------
    spectra
    wave_min
        [um]
    wave_max
        [um]
    area : Quantity
        [m2]
    bandpass : SpectralElement


    Returns
    -------
    counts : u.Quantity array

    """

    if isinstance(wave_min, u.Quantity):
        wave_min = wave_min.to(u.Angstrom).value
    else:
        wave_min *= 1E4

    if isinstance(wave_max, u.Quantity):
        wave_max = wave_max.to(u.Angstrom).value
    else:
        wave_max *= 1E4

    counts = []
    for spec in spectra:
        mask = (spec.model.points[0] > wave_min) * \
               (spec.model.points[0] < wave_max)
        x = spec.model.points[0][mask]
        x = np.append(np.append(wave_min, x), wave_max)
        y = spec.model.lookup_table[mask]
        y = np.append(np.append(spec(wave_min), y), spec(wave_max))

        # flux [ph s-1 cm-2] == y [ph s-1 cm-2 AA-1] * x [AA]
        if isinstance(bandpass, SpectralElement):
            bp = bandpass(x)
            bandpass.model.bounds_error = True
            counts += [np.trapz(bp * y, x)]
        else:
            counts += [np.trapz(y, x)]

    # counts = flux [ph s-1 cm-2]
    counts = 1E4 * np.array(counts)    # to get from cm-2 to m-2
    counts *= u.ph * u.s**-1 * u.m**-2
    if area is not None:
        counts *= utils.quantify(area, u.m ** 2)

    return counts


def new_photons_in_range(spectra, wave_min, wave_max, area, bandpass=None):
    """
    This function intends to supersede the function above. Once tests pass
    It relies in synphot.Observations to return the photon counts

    TODO: Write wrapper functions make_synphot_bandpass and make_synphot_spectra
        to allow a variety of bandpasses and spectra.
    TODO: Delete above function once tests pass


    Parameters
    ----------
    spectra: a synphot spectrum
    wave_min
        [um]
    wave_max
        [um]
    area : Quantity
        [m2]
    bandpass : SpectralElement


    Returns
    -------
    counts : u.Quantity array

    """
    if isinstance(area, u.Quantity):
        area = area.to(u.m**2).value  # if not unit is given, area is assumed in m2

    if isinstance(wave_min, u.Quantity):
        wave_min = wave_min.to(u.Angstrom).value
    else:
        wave_min *= 1E4

    if isinstance(wave_max, u.Quantity):
        wave_max = wave_max.to(u.Angstrom).value
    else:
        wave_max *= 1E4  # if not unit is given, wavelength is assumed in Angstrom

    if isinstance(spectra, list) is False:
        spectra = [spectra]

    if bandpass is None:
        # this makes a bandpass out of wmin and wmax
        mid_point = 0.5*(wave_min + wave_max)
        width = abs(wave_max - wave_min)
        bandpass = SpectralElement(Box1D, amplitude=1, x_0=mid_point, width=width)

    if (bandpass is not None) and (isinstance(bandpass, synphot.spectrum.SpectralElement) is False) :
        # bandpass = make_synphot_bandpass(bandpass) # try to make a synphot bandpass from e.g. filter file
        pass

    counts = []
    for spec in spectra:
        if isinstance(spec, synphot.spectrum.SourceSpectrum) is False:
            #spec = make_synphot_spectrum(spec) # Try to make a synphot spectrum from e.g. file/np.array
            pass

        obs = Observation(spec, bandpass)
        counts.append(obs.countrate(area=area*u.m**2).value)

    counts = np.array(counts) * u.ph * u.s**-1

    return counts



def rebin_spectra(spectra, new_waves):
    """
    Rebin a synphot spectra to a new wavelength grid conserving flux.
    Grid does not need to be linear and can be at higher or lower resolution

    TODO: To resample the spectra at lower resolution a convolution is first needed. Implement!
    TODO: Return the new spectra in the input wavelengths

    Parameters
    ----------
    spectra: a synphot spectra
    new_waves: an array of the output wavelenghts in Angstroms but other units can be
        specified


    Returns
    -------

    A synphot spectra in the new wavelengths

    """
    if isinstance(spectra, synphot.spectrum.SourceSpectrum) is False:
        # spec = make_synphot_spectrum(spec) # Try to make a synphot spectrum from e.g. file/np.array
        raise ValueError(spec, "is not a synphot spectra!")

    if spectra.waveset is None:
        raise ValueError("spectra doesn't have a defined waveset")

    if isinstance(new_waves, u.Quantity):
        new_waves = new_waves.to(u.Angstrom).value

    waves = spectra.waveset.value
    f = np.ones(len(waves))

    filt = SpectralElement(Empirical1D, points=waves, lookup_table=f)
    obs = Observation(spectra, filt, binset=new_waves, force='taper')

    newflux = obs.binflux

    rebin_spec = SourceSpectrum(Empirical1D, points=new_waves,
                                lookup_table=newflux, meta=spectra.meta)

    return rebin_spec


def make_imagehdu_from_table(x, y, flux, pix_scale=1*u.arcsec):

    pix_scale = pix_scale.to(u.deg)
    unit = pix_scale.unit
    x = utils.quantify(x, unit)
    y = utils.quantify(y, unit)

    xpixmin = int(np.floor(np.min(x) / pix_scale))
    ypixmin = int(np.floor(np.min(y) / pix_scale))
    xvalmin = (xpixmin * pix_scale).value
    yvalmin = (ypixmin * pix_scale).value

    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.crpix = [0., 0.]
    the_wcs.wcs.cdelt = [pix_scale.value, pix_scale.value]
    the_wcs.wcs.crval = [xvalmin, yvalmin]
    the_wcs.wcs.cunit = [unit, unit]
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    ypix, xpix = the_wcs.wcs_world2pix(y.to(u.deg), x.to(u.deg), 1)
    yint, xint  = ypix.astype(int), xpix.astype(int)

    image = np.zeros((np.max(xint) + 1, np.max(yint) + 1))
    for ii in range(len(xint)):
        image[xint[ii], yint[ii]] += flux[ii]

    hdu = fits.ImageHDU(data=image)
    hdu.header.extend(the_wcs.to_header())

    return hdu


def scale_imagehdu(imagehdu, waverange, area=None):
    # ..todo: implement this
    # For the moment, all imagehdu must be accompanied by a spectrum in PHOTLAM
    #
    # Future functionality will include scaling here of:
    # ph s-1
    # ph s-1 m-2
    # ph s-1 m-2
    # ph s-1 m-2 um-1
    # ph s-1 m-2 um-1 arcsec-2
    # J s-1 m-2 Hz-1
    # J s-1 m-2 Hz-1 arcsec-2
    # ABMAG
    # ABMAG arcsec-2
    # VEGAMAG
    # VEGAMAG arcsec-2

    if "SPEC_REF" not in imagehdu.header:
        raise ValueError("For this version, an ImageHDU must be associated "
                         "with a spectrum. This will change in the future.")

    return imagehdu

#     unit = extract_unit_from_imagehdu(imagehdu)
#
#     per_unit_area = False
#     if area is None:
#         per_unit_area = True
#         area = 1*u.m**2
#
#     unit, sa_unit = utils.extract_type_from_unit(unit, "solid angle")
#     unit = convert_flux(waverange, 1 * unit, "count", area=area)
#     # [ct] comes out of convert_flux
#     unit *= u.s
#
#     if sa_unit != "":
#         cunit1 = u.deg
#         cunit2 = u.deg
#         if "CUNIT1" in imagehdu.header and "CUNIT2" in imagehdu.header:
#             cunit1 = u.Unit(imagehdu.header["CUNIT1"])
#             cunit2 = u.Unit(imagehdu.header["CUNIT2"])
#         cdelt1 = imagehdu.header["CDELT1"] * cunit1
#         cdelt2 = imagehdu.header["CDELT2"] * cunit2
#
#         pix_angle_area = cdelt1 * cdelt2
#         unit *= (pix_angle_area * sa_unit).si.value
#
#     if per_unit_area is True:
#         unit *= u.m**-2
#
#     zero  = 0 * u.Unit(unit)
#     scale = 1 * u.Unit(unit)
#
#     if "BSCALE" in imagehdu.header:
#         scale *= imagehdu.header["BSCALE"]
#         imagehdu.header["BSCALE"] = 1
#     if "BZERO" in imagehdu.header:
#         zero = imagehdu.header["BZERO"]
#         imagehdu.header["BZERO"] = 0
#
#     imagehdu.data = imagehdu * scale + zero
#     imagehdu.header["BUNIT"] = str(imagehdu.data.unit)
#     imagehdu.header["FLUXUNIT"] = str(imagehdu.data.unit)
#
#     return imagehdu
#
# def extract_unit_from_imagehdu(imagehdu):
#     if "BUNIT" in imagehdu.header:
#         unit = u.Unit(imagehdu.header["BUNIT"])
#     elif "FLUXUNIT" in imagehdu.header:
#         unit = u.Unit(imagehdu.header["FLUXUNIT"])
#     elif isinstance(imagehdu.data, u.Quantity):
#         unit = imagehdu.data.unit
#         imagehdu.data = imagehdu.data.value
#     else:
#         unit = ""
#         warnings.warn("No flux unit found on ImageHDU. Please add BUNIT or "
#                       "FLUXUNIT to the header.")
#
#     return unit


def empty_sky():
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
    from .source import Source
    sky = Source(lam=np.array([0.3, 3.0]), spectra=np.array([0, 0]),
                 x=[0], y=[0], ref=[0], weight=[0])
    return sky
