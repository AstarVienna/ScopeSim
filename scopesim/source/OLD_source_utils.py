import warnings

import numpy as np
from astropy import units as u, constants as c
from astropy.io import ascii as ioascii, fits

from ..OLD_spectral import TransmissionCurve, EmissionCurve
from ..utils import find_file


def _get_stellar_properties(spec_type, cat=None, verbose=False):
    """
    Returns an :class:`astropy.Table` with the list of properties for the
    star(s) in ``spec_type``

    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types
    cat : str, optional
        The filename of a catalogue in a format readable by
        :func:`astropy.io.ascii.read`, e.g. ASCII, CSV. The catalogue should
        contain stellar properties
    verbose : bool
        Print which stellar type is being considered

    Returns
    -------
    props : :class:`astropy.Table` or list of :class:`astropy.Table` py_objects
        with stellar parameters

    """

    if cat is None:
        cat = ioascii.read(find_file("EC_all_stars.csv"))

    if isinstance(spec_type, (list, tuple)):
        return [_get_stellar_properties(i, cat) for i in spec_type]
    else:
        # Check if stellar type is in cat; if not look for the next
        # type in the sequence that is and assign its values
        spt, cls, lum = spec_type[0], int(spec_type[1]), spec_type[2:]
        for _ in range(10):
            if cls > 9:
                cls = 0
                spt = "OBAFGKMLT"["OBAFGKMLT".index(spt)+1]

            startype = spt+str(cls)+lum # was 'star', redefined function star()
            cls += 1

            if startype in cat["Stellar_Type"]:
                break

        else:   # for loop did not find anything
            raise ValueError(spec_type+" doesn't exist in the database")

        n = np.where(cat["Stellar_Type"] == startype.upper())[0][0]
        if verbose:
            print("Returning properties for", startype)

        return cat[n]


def _get_stellar_mass(spec_type):
    """
    Returns a single (or list of) float(s) with the stellar mass(es)

    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types in the normal format: G2V

    Returns
    -------
    mass : float, list
        [Msol]

    """

    props = _get_stellar_properties(spec_type)

    if isinstance(props, (list, tuple)):
        return [prop["Mass"] for prop in props]
    else:
        return props["Mass"]


def _get_stellar_Mv(spec_type):
    """
    Returns a single (or list of) float(s) with the V-band absolute magnitude(s)

    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types

    Returns
    -------
    Mv : float, list

    """

    props = _get_stellar_properties(spec_type)

    if isinstance(props, (list, tuple)):
        return [prop["Mv"] for prop in props]
    else:
        return props["Mv"]


def _get_pickles_curve(spec_type, cat=None, verbose=False):
    """
    Returns the emission curve for a single or list of ``spec_type``, normalised
    to 5556A

    Parameters
    ----------
    spec_type : str, list
        The single (or list) of spectral types (i.e. "A0V" or ["K5III", "B5I"])

    Returns
    -------
    lam : np.array
        a single np.ndarray for the wavelength bins of the spectrum,
    val : np.array (list)
        a (list of) np.ndarray for the emission curve of the spectral type(s)
        relative to the flux at 5556A

    References
    ----------
    Pickles 1998 - DOI: 10.1086/316197

    """
    if cat is None:
        cat = fits.getdata(find_file("EC_pickles.fits"))

    if isinstance(spec_type, (list, tuple)):
        return cat["lam"], [_get_pickles_curve(i, cat)[1] for i in spec_type]
    else:
        # split the spectral type into 3 components and generalise for Pickles
        spt, cls, lum = spec_type[0], int(spec_type[1]), spec_type[2:]
        if lum.upper() == "I":
            lum = "Ia"
        elif lum.upper() == "II":
            lum = "III"
        elif "V" in lum.upper():
            lum = "V"

        for _ in range(10):  # TODO: What does this loop do? (OC)
            if cls > 9:
                cls = 0
                spt = "OBAFGKMLT"["OBAFGKMLT".index(spt)+1]
            startype = spt + str(cls) + lum
            cls += 1

            if startype in cat.columns.names:
                break

        if spec_type != startype and verbose:
            print(spec_type, "isn't in Pickles. Returned", startype)

        try:
            lam, spec = cat["lam"], cat[startype]
        except KeyError:      # Correct? This shouldn't use error handling.
            lam, spec = cat["lam"], cat["M9III"]
        return lam, spec


def _scale_pickles_to_photons(spec_type, mag=0):
    """
    Pull in a spectrum from the Pickles library and scale to V=0 star

    Parameters
    ----------
    spec_type : str, list
        A (list of) spectral type(s), e.g. "A0V" or ["A0V", G2V"]
    mag : float, list, optional
        A (list of) magnitudes for the spectral type(s). Default is 0

    Returns
    -------
    lam, ec : array
        The wavelength bins and the SEDs for the spectral type

    Notes
    -----
    - Vega has a 5556 flux of between 950 and 1000 ph/s/cm2/A. The pickles
    resolution is 5 Ang.
    - Therefore the flux at 5555 should be 5 * 1000 * 10^(-0.4*Mv) ph/s/cm2/bin
    - Pickles catalogue is in units of Flambda [erg/s/cm2/A]
    - Ergo we need to divide the pickels values by lam/0.5556[nm], then rescale
    Regarding the number of photons in the 1 Ang bin at 5556 Ang
    - Bohlin (2014) says F(5556)=3.44x10-9 erg cm-2 s-1 A-1
    - Values range from 3.39 to 3.46 with the majority in range 3.44 to 3.46.
      Bohlin recommends 3.44
    - This results in a photon flux of 962 ph cm-2 s-1 A-1 at 5556 Ang

    """

    if isinstance(spec_type, (list, tuple, np.ndarray)):
        if isinstance(mag, (list, tuple, np.ndarray)):
            if len(mag) != len(spec_type):
                raise ValueError("len(mag) != len(spec_type)")
            mag = list(mag)
        else:
            mag = [mag]*len(spec_type)
    else:
        mag = [mag]

    mag = np.asarray(mag)

    Mv = _get_stellar_Mv(spec_type)
    if not hasattr(Mv, "__len__"):
        Mv = [Mv]

    Mv = np.asarray(Mv)
    lam, ec = _get_pickles_curve(spec_type)
    dlam = (lam[1:] - lam[:-1])
    dlam = np.append(dlam, dlam[-1])

    lam *= 1E-4         # convert to um from Ang

    # Use Bohlin (2014) to determine the photon flux of a mag 0 A0V star
    # at 5556 Ang
    F = 3.44E-9 * u.erg / (u.cm**2 * u.s * u.AA)
    E = c.c*c.h/(5556*u.AA)
    ph0 = (F/E).to(1/(u.s * u.cm**2 * u.AA)).value

    # 5 Ang/bin * ~962 ph/s * (abs mag + apparent mag)

    ph_factor = []
    for i in range(len(mag)):
        tmp = dlam * ph0 * 10**(-0.4*(Mv[i] + mag[i]))
        ph_factor += [tmp]

    # take care of the conversion to ph/s/m2 by multiplying by 1E4
    # TODO: The original type(ec) == (list, tuple) is wrong (should be 'in')
    #   However, correcting it (using idiomatic isinstance) breaks the code!
    #   There must be a bug.
    # Correct code:
    # if isinstance(ec, (list, tuple)):
    #     for i in range(len(ec)):
    if type(ec) == (list, tuple):
        for i in range(len(ec)):
            ec[i] *= (lam/0.5556) * ph_factor[i] * 1E4
    else:
        ec *= (lam/0.5556) * ph_factor[0] * 1E4

    return lam, ec


def BV_to_spec_type(B_V):
    """
    Returns the latest main sequence spectral type(s) for (a) B-V colour

    Parameters
    ----------
    B_V : float, array
        [mag] B-V colour

    Returns
    -------
    spec_types : list
        A list of the spectral types corresponding to the B-V colours

    Examples
    --------
    ::

        >>> BV = np.arange(-0.3, 2.5, 0.5)
        >>> spec_types = BV_to_spec_type(BV)
        >>> print(BV)
        >>> print(spec_types)
        [-0.3  0.2  0.7  1.2  1.7  2.2]
        ['O9V', 'A8V', 'G2V', 'K5V', 'M3V', 'M8V']

    """

    #from scopesim.source import _get_stellar_properties

    spec_type = [spt+str(i)+"V" for spt in "OBAFGKM" for i in range(10)]
    B_V_int = np.array([spt["B-V"] for spt in _get_stellar_properties(spec_type)])

    idx = np.round(np.interp(B_V, B_V_int, np.arange(len(B_V_int)))).astype(int)
    if np.isscalar(idx):
        idx = np.array([idx])
    spec_types = [spec_type[i] for i in idx]

    return spec_types


def mag_to_photons(filter_name, magnitude=0):
    """
    Return the number of photons for a certain filter and magnitude

    Parameters
    ----------
    filter_name : str
        filter name. See scopesim.optics.get_filter_set()
    magnitude : float
        [mag] the source brightness

    Returns
    -------
    flux : float
        [ph/s/m2] Photon flux in the given filter

    See Also
    --------
    :func:`.photons_to_mag`
    :func:`.zero_magnitude_photon_flux`,
    :func:`scopesim.optics.get_filter_set`
    """

    flux_0 = zero_magnitude_photon_flux(filter_name)
    flux = flux_0 * 10**(-0.4 * magnitude)
    return flux


def photons_to_mag(filter_name, photons=1):
    """
    Return the number of photons for a certain filter and magnitude

    Parameters
    ----------
    filter_name : str
        filter name. See scopesim.optics.get_filter_set()
    photons : float
        [ph/s/m2] the integrated photon flux for the filter

    Returns
    -------
    mag : float
        The magnitude of an object with the given photon flux through the filter

    See Also
    --------
    :func:`.photons_to_mag`
    :func:`.zero_magnitude_photon_flux`,
    :func:`scopesim.optics.get_filter_set`

    """

    flux_0 = zero_magnitude_photon_flux(filter_name)
    mag = -2.5 * np.log10(photons / flux_0)
    return mag


def _get_refstar_curve(filename=None, mag=0):
    """
    """
    ## TODO: Can we pre-select a star based on the instrument we're simulating?
    data = ioascii.read(find_file("vega.dat"))
    #data = ioascii.read(find_file("sirius_downsampled.txt"))

    mag_scale_factor = 10**(-mag/2.5)

    ##
    ## this function is expected to return the number of photons of a 0th mag star
    ## for a star brighter than 0th mag, the number of photons needs to be reduced to match a 0th mag star
    lam, spec = data[data.colnames[0]], data[data.colnames[1]]/mag_scale_factor
    return lam, spec


def zero_magnitude_photon_flux(filter_name):
    """
    Return the number of photons for a m=0 star for a certain filter

    Parameters
    ----------
    filter_name : str
        filter name. See scopesim.optics.get_filter_set()

    Notes
    -----
    units in [ph/s/m2]
    """

    if isinstance(filter_name, TransmissionCurve):
        vlam = filter_name.lam
        vval = filter_name.val
    else:
        fname = find_file(filter_name, silent=True)
        if fname is None:
            fname = find_file("TC_filter_" + filter_name + ".dat",
                              silent=True)
            if fname is None:
                raise ValueError("Filter " + filter_name + "cannot be found")

        vraw = ioascii.read(fname)
        vlam = vraw[vraw.colnames[0]]
        vval = vraw[vraw.colnames[1]]

    #lam, vega = _scale_pickles_to_photons("A0V", mag=-0.58)
    ##
    ## we refer here (SimCADO) to the Vega spectrum
    ## (see _get_refstar_curve above).
    lam, vega = _get_refstar_curve(mag=0.)
    filt = np.interp(lam, vlam, vval)

    n_ph = np.sum(vega*filt)

    #print("units in [ph/s/m2]")
    return n_ph


def value_at_lambda(lam_i, lam, val, return_index=False):
    """
    Return the value at a certain wavelength - i.e. val[lam] = x

    Parameters
    ----------
    lam_i : float
        the wavelength of interest
    lam : np.ndarray
        an array of wavelengths
    val : np.ndarray
        an array of values
    return_index : bool, optional
        If True, the index of the wavelength of interest is returned
        Default is False
    """

    i0 = np.where(lam <= lam_i)[0][-1]
    i1 = np.where(lam > lam_i)[0][0]

    lam_x = np.array([lam[i0], lam_i, lam[i1]])
    val_i = np.interp(lam_x, lam, val)

    if return_index:
        return i0
    else:
        return val_i[1]


def scale_spectrum(lam, spec, mag, filter_name="Ks", return_ec=False):
    """
    Scale a spectrum to be a certain magnitude

    Parameters
    ----------
    lam : np.ndarray
        [um] The wavelength bins for spectrum
    spec : np.ndarray
        The spectrum to be scaled into [ph/s/m2] for the given broadband filter
    mag : float
        magnitude of the source
    filter_name : str, TransmissionCurve, optional
           Any filter name from SimCADO or a
           :class:`~.scopesim.spectral.TransmissionCurve` object
           (see :func:`~.scopesim.optics.get_filter_set`)
    return_ec : bool, optional
        If True, a :class:`scopesim.spectral.EmissionCurve` object is returned.
        Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.ndarray
        [ph/s/m2] The spectrum scaled to the specified magnitude

    If return_ec == True, a :class:`scopesim.spectral.EmissionCurve` is returned

    See Also
    --------
    :class:`scopesim.spectral.TransmissionCurve`,
    :func:`scopesim.optics.get_filter_curve`,
    :func:`scopesim.optics.get_filter_set`,
    :func:`scopesim.source.SED`,
    :func:`scopesim.source.stars`

    Examples
    --------

    Scale the spectrum of a G2V star to J=25::

        >>> lam, spec = scopesim.source.SED("G2V")
        >>> lam, spec = scopesim.source.scale_spectrum(lam, spec, 25, "J")

    Scale the spectra for many stars to different H-band magnitudes::

        >>> from scopesim.source.source_utils import scale_spectrum
        >>> from scopesim.source.templates import SED
        >>>
        >>> star_list = ["A0V", "G2V", "M5V", "B6III", "O9I", "M2IV"]
        >>> magnitudes = [ 20,  25.5,  29.1,      17,  14.3,   22   ]
        >>> lam, spec = SED(star_list)
        >>> lam, spec = scale_spectrum(lam, spec, magnitudes, "H")

    Re-scale the above spectra to the same magnitudes in Pa-Beta::

        >>> # Find which filters are in the scopesim/data directory
        >>>
        >>> import scopesim.optics as sim_op
        >>> print(sim_op.get_filter_set())
        ['B', 'BrGamma', 'CH4_169', 'CH4_227', 'FeII_166', 'H', 'H2O_204',
            'H2_212', 'Hcont_158', 'I', 'J', 'K', 'Ks', 'NH3_153', 'PaBeta',
            'R', 'U', 'V', 'Y', 'z']
        >>>
        >>> lam, spec = scale_spectrum(lam, spec, magnitudes, "PaBeta")

    """
    # The following was part of docstring. The example does not work, because
    # the new filter is not calibrated.
    #
    # Create a tophat filter and rescale to magnitudes in that band:
    #
    #     >>> # first make a transmission curve for the filter
    #     >>>
    #     >>> from scopesim.spectral import TransmissionCurve
    #     >>> filt_lam   = np.array([0.3, 1.09, 1.1, 1.15, 1.16, 3.])
    #     >>> filt_trans = np.array([0.,  0.,   1.,  1.,   0.,   0.])
    #     >>> new_filt   = TransmissionCurve(lam=filt_lam, val=filt_trans)
    #     >>>
    #     >>> lam, spec = scale_spectrum(lam, spec, magnitudes, new_filt)

    from ..optics.OLD_optics_utils import get_filter_curve

    mag = np.asarray(mag)

    # Number of photons corresponding to desired apparent magnitude mag
    ideal_phs = zero_magnitude_photon_flux(filter_name) * 10**(-0.4 * mag)
    if isinstance(ideal_phs, (int, float)):
        ideal_phs = [ideal_phs]

    if len(np.shape(spec)) == 1:
        spec = [spec]

    # Convert spectra to EmissionCurves
    curves = [EmissionCurve(lam=lam, val=sp, area=1, units="ph/s/m2")
              for sp in spec]

    if isinstance(filter_name, TransmissionCurve):
        filt = filter_name
    else:
        fname = find_file(filter_name)
        if fname is not None:
            filt = TransmissionCurve(filename=fname)
        else:
            filt = get_filter_curve(filter_name)

    # Rescale the spectra
    for i in range(len(curves)):
        tmp = curves[i] * filt
        obs_ph = tmp.photons_in_range()
        scale_factor = ideal_phs[i] / obs_ph
        curves[i] *= scale_factor

    # Return in desired format
    if return_ec:
        if len(curves) > 1:
            return curves
        else:
            return curves[0]
    else:
        if len(curves) > 1:
            return curves[0].lam, np.array([curve.val for curve in curves])
        else:
            return curves[0].lam, curves[0].val


def scale_spectrum_sb(lam, spec, mag_per_arcsec, pix_res=0.004,
                      filter_name="Ks", return_ec=False):
    """
    Scale a spectrum to be a certain magnitude per arcsec2

    Parameters
    ----------
    lam : np.ndarray
        [um] The wavelength bins for spectrum
    spec : np.ndarray
        The spectrum to be scaled into [ph/s/m2] for the given broadband filter
    mag_per_arcsec : float
        [mag/arcsec2] surface brightness of the source
    pix_res : float
        [arcsec] the pixel resolution
    filter_name : str, TransmissionCurve
        Any filter name from SimCADO or a
        :class:`~.scopesim.spectral.TransmissionCurve` object
        (see :func:`~.scopesim.optics.get_filter_set`)
    return_ec : bool, optional
        If True, a :class:`scopesim.spectral.EmissionCurve` object is returned.
        Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2/pixel] The spectrum scaled to the specified magnitude

    See Also
    --------

    """

    curve = scale_spectrum(lam, spec, mag_per_arcsec, filter_name,
                           return_ec=True)
    curve.val *= pix_res**2
    curve.params["pix_res"] = pix_res

    if return_ec:
        return curve
    else:
        return curve.lam, curve.val


def get_lum_class_params(lum_class="V", cat=None):
    """
    Returns a table with parameters for a certain luminosity class

    Parameters
    ----------
    lum_class : str, optional
        Default is the main sequence ("V")

    Returns
    -------
    :class:`astropy.table.Table` object

    """

    import astropy.table as tbl

    if cat is None:
        cat = ioascii.read(find_file("EC_all_stars.csv"))

    t = []
    for row in cat:
        spt = row["Stellar_Type"]
        if spt[0] in "OBAFGKM" and \
           spt[-len(lum_class):] == lum_class and \
           len(spt) == 2 + len(lum_class):
            t += [row.data]

    t = tbl.Table(data=np.array(t), names=cat.colnames)

    return t


def get_nearest_spec_type(value, param="B-V", cat=None):
    """
    Return the spectral type of the star with the closest parameter value

    Compares values given for a certain stellar parameter and returns the
    spectral type which matches the best. In case several spectral types have
    the same value, only the first spectral type is returned

    Acceptable parameters are:
    "Mass" : [Msun]
    "Luminosity" : [Lsun]
    "Radius" : [Rsun]
    "Temp" : [K]
    "B-V" : [mag]
    "Mv" : [mag]
    "BC(Temp)" : [Corr]
    "Mbol" : [mag]

    Parameters
    ----------
    value : float, array
        The value that the spectral type should have
    param : str, optional
        Default is "B-V". The column to be searched.
    cat : astropy.Table, optional
        The catalogue to use. Default is in the scopesim/data directory

    Returns
    -------
    spec_type : str, list
        a value/list of strings corresponding to the spectral types which best
        fit to the given values

    """

    if cat is None:
        cat = ioascii.read(find_file("EC_all_stars.csv"))

    if isinstance(value, (np.ndarray, list, tuple)):
        spt = []
        for val in value:
            spt += [get_nearest_spec_type(val, param, cat)]

        return spt

    col = cat[param]
    i = np.argmin(np.abs(col-value))
    spec_type = cat["Stellar_Type"][i]

    return spec_type


def spectrum_sum_over_range(lam, flux, lam_min=None, lam_max=None):
    """
    Sum spectrum over range lam_min to lam_max

    Parameters
    ----------
    lam : float, array
        wavelength array of spectrum

    flux : float, array
        flux array of spectrum [ph/s/m2/bin]

    lam_min, lam_max : float
        wavelength limits of range over which the spectrum is summed. If None,
        the spectrum is summed over its definition range


    Returns
    -------
    spec_photons : float
        number of photons within lam_min and lam_max [ph/s/m2]

    """

    if lam_min is None:
        lam_min = lam[0]
    if lam_max is None:
        lam_max = lam[-1]

    if lam_max < lam_min:
        raise ValueError("lam_max < lam_min")

    # Check if the slice limits are within the spectrum wavelength range
    dlam = lam[1] - lam[0]
    if (lam_min > lam[-1] + dlam/2) or (lam_max < lam[0] - dlam/2):
        print((lam_min, lam_max), (lam[0], lam[-1]))
        warnings.warn("lam_min or lam_max outside wavelength range" +
                      " of spectra. Returning 0 photons for this range")
        return np.array([0])

    # find the closest indices imin, imax that match the limits
    imin = np.argmin(np.abs(lam - lam_min))
    imax = np.argmin(np.abs(lam - lam_max))

    # Treat edge bins: Since lam[imin] < lam_min < lam_max < lam[imax], we have
    # to subtract part of the outer bins
    dlam = lam[1] - lam[0]
    if np.size(np.shape(flux)) == 1: # 1D case
        spec_photons = np.sum(flux[imin:(imax + 1)]) \
                        - flux[imin] * (0.5 + (lam_min - lam[imin])/dlam) \
                        - flux[imax] * (0.5 - (lam_max - lam[imax])/dlam)

    elif np.size(np.shape(flux)) == 2: # 2D case
        spec_photons = np.sum(flux[:, imin:(imax + 1)], axis=1) \
                        - flux[:, imin] * (0.5 + (lam_min - lam[imin])/dlam) \
                        - flux[:, imax] * (0.5 - (lam_max - lam[imax])/dlam)
    else:
        spec_photons = 0

    return spec_photons


def _rebin(img, bpix):
    """Rebin image img by block averaging bpix x bpix pixels"""

    xedge = np.shape(img)[0] % bpix
    yedge = np.shape(img)[1] % bpix
    img_block = img[xedge:, yedge:]

    binim = np.reshape(img_block,
                       (int(img_block.shape[0]/bpix), bpix,
                        int(img_block.shape[1]/bpix), bpix))
    binim = np.mean(binim, axis=3)
    binim = np.mean(binim, axis=1)
    return binim