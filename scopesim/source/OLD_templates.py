import os
from glob import glob

import numpy as np
from astropy import units as u
from astropy.io import ascii as ioascii, fits
from scipy.ndimage import interpolation as spi

from .. import rc, utils
from ..OLD_spectral import EmissionCurve
from ..utils import find_file
from .OLD_source import Source
from .OLD_source_utils import scale_spectrum, _scale_pickles_to_photons, \
    _get_stellar_mass, _get_stellar_Mv, scale_spectrum_sb, _rebin


def get_SED_names(path=None):
    """
    Return a list of the SEDs installed in the package directory

    Looks for files that follow the naming convention ``SED_<name>.dat``.
    For example, SimCADO contains an SED for an elliptical galaxy named
    ``SED_elliptical.dat``

    Parameters
    ----------
    path : str, optional
        Directory to look in for filters

    Returns
    -------
    sed_names : list
        A list of names for the SED files available

    Examples
    --------
    Names returned here can be used with the function :func:`.SED` to call up
    ::

        >>> from scopesim.source import SED, get_SED_names
        >>> print(get_SED_names())
        ['elliptical', 'interacting', 'spiral', 'starburst', 'ulirg']
        >>> SED("spiral")
        (array([ 0.3  ,  0.301,  0.302, ...,  2.997,  2.998,  2.999]),
         array([        0.        ,         0.        ,  26055075.98709349, ...,
                  5007498.76444208,   5000699.21993188,   4993899.67542169]))

    See Also
    --------
    :func:`.SED`

    """
    if path is None:
        path = rc.__data_dir__
    sed_names = [i.replace(".dat", "").split("SED_")[-1] \
                                for i in glob(os.path.join(path, "SED_*.dat"))]

    sed_names += ["All stellar spectral types (e.g. G2V, K0III)"]
    return sed_names


def SED(spec_type, filter_name="V", magnitude=0.):
    """
    Return a scaled SED for a star or type of galaxy

    The SED can be for stellar spectra of galacty spectra. It is best not to mix
    the two types when calling ``SED()``. Either provide a list of stellar types,
    e.g. ["G2V", "A0V"], of a list of galaxy types, e.g. ["elliptical", "starburst"]

    To get the list of galaxy types that are installed, call get_SED_names().
    All stellar types from the Pickles (1998) catalogue are available.

    Parameters
    ----------
    spec_type : str, list
        The spectral type of the star(s) - from the Pickles 1998 catalogue
        The names of a galaxy spectrum - see get_SED_names()
    filter_name : str, optional
        Default is "V". Any filter in the scopesim/data directory can be used,
        or the user can specify a file path to an ASCII file for the filter
    magnitude : float, list, optional
        Apparent magnitude of the star. Default is 0.

    Returns
    -------
    lam : np.ndarray
        [um] The centre of each 5 Ang bin along the spectral axis
    val : np.ndarray
        [ph/s/m2/bin] The photon flux of the star in each bin


    Examples
    --------

    Get the SED and the wavelength bins for a J=0 A0V star

        >>> from scopesim.source import SED
        >>> lam, spec = SED("A0V", "J", 0)

    Get the SED for a generic starburst galaxy

        >>> lam, spec = SED("starburst")

    Get the SEDs for several spectral types with different magnitudes

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from scopesim.source import SED

        lam, spec = SED(spec_type=["A0V", "G2V"],
                            filter_name="PaBeta",
                            magnitude=[15, 20])

        plt.plot(lam, spec[0], "blue", label="Vega")
        plt.plot(lam, spec[1], "orange", label="G2V")
        plt.semilogy(); plt.legend(); plt.show()



    Notes
    -----
    Original flux units for the stellar spectra are in [ph/s/m2/AA], so we
    multiply the flux by 5 to get [ph/s/m2/bin]. Therefore divide by 5*1E4 if
    you need the flux in [ph/s/cm2/Angstrom]

    """

    if isinstance(spec_type, (tuple, list, np.ndarray)):
        spec_type = list(spec_type)
        if np.isscalar(magnitude):
            magnitude = [magnitude]*len(spec_type)
    elif isinstance(spec_type, str):
        spec_type = [spec_type]

    if isinstance(magnitude, (list, tuple)):
        magnitude = np.asarray(magnitude)

    # Check if any of the names given are in the package directory
    gal_seds = get_SED_names()
    if np.any([i in gal_seds for i in spec_type]):
        galflux = []
        for gal in spec_type:
            data = ioascii.read(find_file("/data/SED_"+gal+".dat"))
            galflux += [data[data.colnames[1]]]
            galflux = np.asarray(galflux)
        lam = data[data.colnames[0]]

        lam, galflux = scale_spectrum(lam=lam, spec=galflux, mag=magnitude,
                                      filter_name=filter_name)
        return lam, galflux

    else:
        lam, starflux = _scale_pickles_to_photons(spec_type)
        lam, starflux = scale_spectrum(lam=lam, spec=starflux, mag=magnitude,
                                       filter_name=filter_name)

        return lam, starflux


def empty_sky():
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
    sky = Source(lam=np.linspace(0.3, 3.0, 271),
                 spectra=np.zeros((1, 271)),
                 x=[0], y=[0], ref=[0], weight=[0])
    return sky


def star_grid(n, mag_min, mag_max, filter_name="Ks", separation=1,
              spec_type="A0V"):
    """
    Creates a square grid of A0V stars at equal magnitude intervals

    Parameters
    ----------
    n : float
        the number of stars in the grid
    mag_min, mag_max : float
        [vega mag] the minimum (brightest) and maximum (faintest) magnitudes for
        stars in the grid
    filter_name : str
        any filter that is in the SimCADO package directory.
        See ``scopesim.optics.get_filter_set()``
    separation : float, optional
        [arcsec] an average speration between the stars in the grid can be
        specified. Default is 1 arcsec
    spec_type : str, optional
        the spectral type of the star, e.g. "A0V", "G5III"

    Returns
    -------
    source : ``scopesim.Source``

    Notes
    -----
    The units of the A0V spectrum in ``source`` are [ph/s/m2/bin].
    The weight values are the scaling factors to bring a V=0 A0V spectrum down
    to the required magnitude for each star.

    """

    if isinstance(mag_min, (list, tuple, np.ndarray)):
        mags = np.asarray(mag_min)
    else:
        if mag_min < mag_max:
            mags = np.linspace(mag_min, mag_max, n)
        elif mag_min > mag_max:
            mags = np.linspace(mag_max, mag_min, n)
        elif mag_min == mag_max:
            mags = np.ones(n) * mag_min

    side_len = int(np.sqrt(n)) + (np.sqrt(n) % 1 > 0)

    x = separation * (np.arange(n) % side_len - (side_len - 1) / 2)
    y = separation * (np.arange(n)// side_len - (side_len - 1) / 2)

    lam, spec = SED(spec_type, filter_name=filter_name, magnitude=0)
    if isinstance(spec_type, (list, tuple)):
        ref = np.arange(len(spec_type))
    else:
        ref = np.zeros((n))
    weight = 10**(-0.4*mags)

    units = "ph/s/m2"

    src = Source(lam=lam, spectra=spec,
                 x=x, y=y,
                 ref=ref, weight=weight,
                 units=units)

    return src


def star(spec_type="A0V", mag=0, filter_name="Ks", x=0, y=0, **kwargs):
    """
    Creates a scopesim.Source object for a star with a given magnitude

    This is just the single star variant for ``scopesim.source.stars()``

    Parameters
    ----------
    spec_type : str
        the spectral type of the star, e.g. "A0V", "G5III"
    mag : float
        magnitude of star
    filter_name : str
        Filter in which the magnitude is given. Can be the name of any filter
        curve file in the scopesim/data folder, or a path to a custom ASCII file
    x, y : float, int, optional
        [arcsec] the x,y position of the star on the focal plane


    Keyword arguments
    -----------------
    Passed to the ``scopesim.Source`` object. See the docstring for this object.

    pix_unit : str
        Default is "arcsec". Acceptable are "arcsec", "arcmin", "deg", "pixel"
    pix_res : float
        [arcsec] The pixel resolution of the detector. Useful for surface
        brightness calculations

    Returns
    -------
    source : ``scopesim.Source``

    See Also
    --------
    .stars()

    """

    thestar = stars([spec_type], [mag], filter_name, [x], [y], **kwargs)
    return thestar


def stars(spec_types=("A0V"), mags=(0), filter_name="Ks",
          x=None, y=None, **kwargs):
    """
    Creates a scopesim.Source object for a bunch of stars.

    Parameters
    ----------
    spec_types : str, list of strings
        the spectral type(s) of the stars, e.g. "A0V", "G5III"
        Default is "A0V"
    mags : float, array
        [mag] magnitudes of the stars.
    filter_name : str,
        Filter in which the magnitude is given. Can be the name of any filter
        curve file in the scopesim/data folder, or a path to a custom ASCII file
    x, y : arrays
        [arcsec] x and y coordinates of the stars on the focal plane


    Keyword arguments
    -----------------
    Passed to the ``scopesim.Source`` object. See the docstring for this object.

    pix_unit : str
        Default is "arcsec". Acceptable are "arcsec", "arcmin", "deg", "pixel"
    pix_res : float
        [arcsec] The pixel resolution of the detector. Useful for surface
        brightness calculations

    Returns
    -------
    source : ``scopesim.Source``


    Examples
    --------

    Create a ``Source`` object for a random group of stars

        >>> import numpy as np
        >>> from scopesim.source import stars
        >>>
        >>> spec_types = ["A0V", "G2V", "K0III", "M5III", "O8I"]
        >>> ids = np.random.randint(0,5, size=100)
        >>> star_list = [spec_types[i] for i in ids]
        >>> mags = np.random.normal(20, 3, size=100)
        >>>
        >>> src = stars(spec_types, mags, filter_name="Ks")

    If we don't specify any coordinates all stars have the position (0, 0).
    **All positions are in arcsec.**
    There are two possible ways to add positions. If we know them to begin with
    we can add them when generating the source full of stars

        >>> x, y = np.random.random(-20, 20, size=(100,2)).tolist()
        >>> src = stars(star_list, mags, filter_name="Ks", x=x, y=y)

    Or we can add them to the ``Source`` object directly (although, there are
    less checks to make sure the dimensions match here):

        >>> src.x, src.y = x, y


    """

    if isinstance(spec_types, str):
        spec_types = [spec_types]

    if isinstance(mags, (int, float)):
        mags = [mags] * len(spec_types)

    if len(mags) > 1  and len(spec_types) == 1:
        spec_types *= len(mags)
    elif len(mags) != len(spec_types):
        raise ValueError("len(mags) != len(spec_types)")

    mags = np.array(mags)

    if x is None:
        x = np.zeros(len(mags))
    if y is None:
        y = np.zeros(len(mags))

    # only pull in the spectra for unique spectral types

    # assign absolute magnitudes to stellar types in cluster
    unique_types = np.unique(spec_types)
    lam, spec = SED(unique_types, filter_name=filter_name,
                    magnitude=[0]*len(unique_types))

    # get the references to the unique stellar types
    ref_dict = {i : j for i, j in zip(unique_types,
                                      np.arange(len(unique_types)))}
    if isinstance(spec_types, (list, tuple, np.ndarray)):
        ref = np.array([ref_dict[i] for i in spec_types])
    else:
        ref = np.zeros(len(mags))

    weight = 10**(-0.4*mags)

    units = "ph/s/m2"

    src = Source(lam=lam, spectra=spec,
                 x=x, y=y,
                 ref=ref, weight=weight,
                 units=units, **kwargs)

    src.info["object"] = "stars"
    src.info["spec_types"] = spec_types
    src.info["magnitudes"] = mags
    src.info["filter_name"] = filter_name

    return src


def source_1E4_Msun_cluster(distance=50000, half_light_radius=1):
    """
    Generate a source object for a 10^4 solar mass cluster

    Parameters
    ----------
    distance : float
        [pc] distance to the cluster
    half_light_radius : float
        [pc] half light radius of the cluster
    mass : float
        [Msun] If you'd like a different size cluster

    Returns
    -------
    src : scopesim.Source

    See Also
    --------
    .cluster()

    """
    # IMF is a realisation of stellar masses drawn from an initial mass
    # function (TODO: which one?) summing to 1e4 M_sol.
    fname = find_file("IMF_1E4.dat")
    imf = np.loadtxt(fname)

    # Assign stellar types to the masses in imf using list of average
    # main-sequence star masses:
    stel_type = [i + str(j) + "V" for i in "OBAFGKM" for j in range(10)]
    mass = _get_stellar_mass(stel_type)
    ref = utils.nearest(mass, imf)
    thestars = [stel_type[i] for i in ref] # was stars, redefined function name

    # assign absolute magnitudes to stellar types in cluster
    unique_ref = np.unique(ref)
    unique_type = [stel_type[i] for i in unique_ref]
    unique_Mv = _get_stellar_Mv(unique_type)

    # Mv_dict = {i : float(str(j)[:6]) for i, j in zip(unique_type, unique_Mv)}
    ref_dict = {i: j for i, j in zip(unique_type, np.arange(len(unique_type)))}

    # find spectra for the stellar types in cluster
    lam, spectra = _scale_pickles_to_photons(unique_type)

    # this one connects the stars to one of the unique spectra
    stars_spec_ref = [ref_dict[i] for i in thestars]

    # absolute mag + distance modulus
    m = np.array([unique_Mv[i] for i in stars_spec_ref])
    m += 5 * np.log10(distance) - 5

    # set the weighting
    weight = 10**(-0.4*m)

    # draw positions of stars: cluster has Gaussian profile
    distance *= u.pc
    half_light_radius *= u.pc
    hwhm = (half_light_radius/distance*u.rad).to(u.arcsec).value
    sig = hwhm / np.sqrt(2 * np.log(2))

    x = np.random.normal(0, sig, len(imf))
    y = np.random.normal(0, sig, len(imf))

    src = Source(lam=lam, spectra=spectra, x=x, y=y, ref=stars_spec_ref,
                 weight=weight, units="ph/s/m2")

    return src


def cluster(mass=1E3, distance=50000, half_light_radius=1):
    """
    Generate a source object for a cluster

    The cluster distribution follows a gaussian profile with the
    ``half_light_radius`` corresponding to the HWHM of the distribution. The
    choice of stars follows a Kroupa IMF, with no evolved stars in the mix. Ergo
    this is more suitable for a young cluster than an evolved custer

    Parameters
    ----------
    mass : float
        [Msun] Mass of the cluster (not number of stars). Max = 1E5 Msun
    distance : float
        [pc] distance to the cluster
    half_light_radius : float
        [pc] half light radius of the cluster

    Returns
    -------
    src : scopesim.Source

    Examples
    --------

    Create a ``Source`` object for a young open cluster with half light radius
    of around 0.2 pc at the galactic centre and 100 solar masses worth of stars:

        >>> from scopesim.source import cluster
        >>> src = cluster(mass=100, distance=8500, half_light_radius=0.2)


    """
    # IMF is a realisation of stellar masses drawn from an initial mass
    # function (TODO: which one?) summing to 1e4 M_sol.
    if mass <= 1E4:
        fname = find_file("IMF_1E4.dat")
        imf = np.loadtxt(fname)
        imf = imf[0:int(mass/1E4 * len(imf))]
    elif mass > 1E4 and mass < 1E5:
        fname = find_file("IMF_1E5.dat")
        imf = np.loadtxt(fname)
        imf = imf[0:int(mass/1E5 * len(imf))]
    else:
        raise ValueError("Mass too high. Must be <10^5 Msun")

    # Assign stellar types to the masses in imf using list of average
    # main-sequence star masses:
    stel_type = [i + str(j) + "V" for i in "OBAFGKM" for j in range(10)]
    masses = _get_stellar_mass(stel_type)
    ref = utils.nearest(masses, imf)
    thestars = [stel_type[i] for i in ref] # was stars, redefined function name

    # assign absolute magnitudes to stellar types in cluster
    unique_ref = np.unique(ref)
    unique_type = [stel_type[i] for i in unique_ref]
    unique_Mv = _get_stellar_Mv(unique_type)

    # Mv_dict = {i : float(str(j)[:6]) for i, j in zip(unique_type, unique_Mv)}
    ref_dict = {i : j for i, j in zip(unique_type, np.arange(len(unique_type)))}

    # find spectra for the stellar types in cluster
    lam, spectra = _scale_pickles_to_photons(unique_type)

    # this one connects the stars to one of the unique spectra
    stars_spec_ref = [ref_dict[i] for i in thestars]

    # absolute mag + distance modulus
    m = np.array([unique_Mv[i] for i in stars_spec_ref])
    m += 5 * np.log10(distance) - 5

    # set the weighting
    weight = 10**(-0.4*m)

    # draw positions of stars: cluster has Gaussian profile
    distance *= u.pc
    half_light_radius *= u.pc
    hwhm = (half_light_radius/distance*u.rad).to(u.arcsec).value
    sig = hwhm / np.sqrt(2 * np.log(2))

    x = np.random.normal(0, sig, len(imf))
    y = np.random.normal(0, sig, len(imf))

    src = Source(lam=lam, spectra=spectra, x=x, y=y, ref=stars_spec_ref,
                 weight=weight, units="ph/s/m2")

    src.info["object"] = "cluster"
    src.info["total_mass"] = mass
    src.info["masses"] = imf
    src.info["half_light_radius"] = half_light_radius
    src.info["hwhm"] = hwhm
    src.info["distance"] = distance
    src.info["stel_type"] = stel_type

    return src


def source_from_image(images, lam, spectra, plate_scale, oversample=1,
                      units="ph/s/m2", flux_threshold=0,
                      center_offset=(0, 0),
                      conserve_flux=True,
                      **kwargs):
    """
    Create a Source object from an image or a list of images.

    .. note::
        ``plate_scale`` is the original plate scale of the images. If this is
        not the same as the plate scale of the ``Detector``
        then you will need to specify oversample to interpolate between the two
        scales. I.e.  oversample = Image plate scale / Detector plate scale


    Parameters
    ----------
    images : np.ndarray, list
        A single or list of np.ndarrays describing where the flux is coming from
        The spectrum for each pixel in the image is weighted by the pixel value.

    lam : np.ndarray
        An array contains the centres of the wavelength bins for the spectra

    spectra : np.ndarray
        A (n,m) array with n spectra, each with m bins

    plate_scale : float
        [arcsec] The plate scale of the images in arcseconds (e.g. 0.004"/pixel)

    oversample : int
        The factor with which to oversample the image. Each image pixel is split
        into (oversample)^2 individual point files.

    units : str, optional
        The energy units of the spectra. Default is [ph/s/m2]

    flux_threshold : float, optional
        If there is noise in the image, set threshold to the noise limit so that
        only real photon files are extracted. Default is 0.

    center_offset : (float, float)
        [arcsec] If the centre of the image is offset, add this offset to (x,y)
        coordinates.

    conserve_flux : bool, optional
        If True, when the image is rescaled, flux is conserved
        i.e. np.sum(image) remains constant
        If False, the maximum value of the image stays constant after rescaling
        i.e. np.max(image) remains constant


    Keyword arguments
    -----------------
    Passed to the ``scopesim.Source`` object. See the docstring for this object.

    pix_unit : str
        Default is "arcsec". Acceptable are "arcsec", "arcmin", "deg", "pixel"

    pix_res : float
        [arcsec] The pixel resolution of the detector. Useful for surface
        brightness calculations


    Returns
    -------
    src : :class:`.Source` object


    Examples
    --------
    To create a :class:`.Source` object we need an image that describes the
    spatial distribution of the object of interest and spectrum. For the sake of
    ease we will assign a generic elliptical galaxy spectrum to the image.::

        >>> from astropy.io import fits
        >>> from scopesim.source import SED, source_from_image
        >>>
        >>> im = fits.getdata("galaxy.fits")
        >>> lam, spec = SED("elliptical")
        >>> src = source_from_image(im, lam, spec,
                                    plate_scale=0.004)

    **Note** Here we have assumed that the plate scale of the image is the same
    as the MICADO wide-field mode, i.e. 0.004 arcseconds. If the image is from a
    real observation, or it was generated with a different pixel scale, we will
    need to tell SimCADO about this::

        >>> src = source_from_image(im, lam, spec,
                                    plate_scale=0.01,
                                    oversample=2.5)

    If the image is from real observations, chances are good that the background
    flux is higher than zero. We can set a ``threshold`` in order to tell
    SimCADO to ignore all pixel with values below the background level::

        >>> src = source_from_image(im, lam, spec,
                                    plate_scale=0.01,
                                    oversample=2.5,
                                    flux_threshold=0.2)

    Finally, if the image centre is not the centre of the observation, we can
    shift the image relative to the MICADO field of view. The units for the
    offset are [arcsec]::

        >>> src = source_from_image(im, lam, spec,
                                    plate_scale=0.01,
                                    oversample=2.5,
                                    flux_threshold=0.2,
                                    center_offset=(10,-15))

    """

    if isinstance(images, (list, tuple)):
        srclist = [source_from_image(images[i], lam, spectra[i, :], plate_scale,
                                     oversample, units, flux_threshold,
                                     center_offset)
                   for i in range(len(images))]
        fullsrc = srclist[0]
        for src in srclist[1:]:
            fullsrc += src
        return fullsrc

    else:
        # if not isinstance(oversample, int):
        #    raise ValueError("Oversample must be of type 'int'")

        if isinstance(images, str) and images.split(".")[-1].lower() == "fits":
            images = fits.getdata(images)

        # im = images
        # y_cen, x_cen = np.array(im.shape) / 2 + np.array(center_offset)
        # # x_cen, y_cen = np.array(im.shape) / 2 + np.array(center_offset)
        # # x_i, y_i = np.where(im > flux_threshold)
        # y_i, x_i = np.where(im > flux_threshold)

        # x = (x_i - x_cen) * plate_scale
        # y = (y_i - y_cen) * plate_scale
        # # weight = im[x_i, y_i]
        # weight = im[y_i, x_i]

        # i = oversample
        # oset = np.linspace(-0.5, 0.5, 2*i+1)[1:2*i:2] * plate_scale

        # x_list, y_list, w_list = [], [], []
        # for i in oset:
            # for j in oset:
                # x_list += (x + i).tolist()
                # y_list += (y + j).tolist()
                # w_list += (weight / oversample**2).tolist()
        # x, y, weight = np.array(x_list), np.array(y_list), np.array(w_list)

        if oversample != 1:
            img = spi.zoom(images, oversample, order=3).astype(np.float32)
            scale_factor = np.sum(images)/np.sum(img)
            if conserve_flux:
                img *= scale_factor
        else:
            img = images
            scale_factor = 1

        # Ugly stripes are fixed - KL - 22.08.2017
        y_cen, x_cen = np.array(img.shape) // 2 + 0.5
        #y_cen, x_cen = np.array(img.shape) / 2
        y_i, x_i = np.where(img > flux_threshold * scale_factor)

        pix_res = plate_scale / oversample
        x = (x_i - x_cen) * pix_res + center_offset[0]
        y = (y_i - y_cen) * pix_res + center_offset[1]

        weight = img[y_i, x_i]
        ref = np.zeros(len(x))

        src = Source(lam=lam, spectra=spectra, x=x, y=y, ref=ref, weight=weight,
                     units=units, **kwargs)

        return src


def flat_spectrum(mag, filter_name="Ks", return_ec=False):
    """
    Return a flat spectrum scaled to a certain magnitude

    Parameters
    ----------
    mag : float
        [mag] magnitude of the source
    filter_name : str, TransmissionCurve, optional
        str - filter name. See ``scopesim.optics.get_filter_set()``. Default: "Ks"
        TransmissionCurve - output of ``scopesim.optics.get_filter_curve()``
    return_ec : bool, optional
        If True, a scopesim.spectral.EmissionCurve object is returned.
        Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2/arcsec] The spectrum scaled to the specified magnitude

    """
    lam = np.arange(0.3, 3.0, 0.01)
    spec = np.ones(len(lam))

    # .. todo: mag_per_arcsec undefined? (OC)
    if return_ec:
        curve = scale_spectrum(lam, spec, mag, filter_name,
                               return_ec)
        return curve
    else:
        lam, spec = scale_spectrum(lam, spec, mag, filter_name,
                                   return_ec)
        return lam, spec


def flat_spectrum_sb(mag_per_arcsec, filter_name="Ks", pix_res=0.004,
                     return_ec=False):
    """
    Return a flat spectrum for a certain magnitude per arcsec

    Parameters
    ----------
    mag_per_arcsec : float
        [mag/arcsec2] surface brightness of the source
    filter_name : str, TransmissionCurve, optional
        str - filter name. See ``scopesim.optics.get_filter_set()``. Default: "Ks"
        TransmissionCurve - output of ``scopesim.optics.get_filter_curve()``
    pix_res : float
        [arcsec] the pixel resolution. Default is 4mas (i.e. 0.004)
    return_ec : bool, optional
        Default is False. If True, a scopesim.spectral.EmissionCurve object is
        returned.

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2/arcsec] The spectrum scaled to the specified magnitude

    """

    lam = np.arange(0.3, 3.0, 0.01)
    spec = np.ones(len(lam))

    if return_ec:
        curve = scale_spectrum_sb(lam, spec, mag_per_arcsec, pix_res,
                                  filter_name, return_ec)
        return curve
    else:
        lam, spec = scale_spectrum_sb(lam, spec, mag_per_arcsec, pix_res,
                                      filter_name, return_ec)
        return lam, spec


def sie_grad(x, y, par):
    """
    Compute the deflection of an SIE (singular isothermal ellipsoid) potential

    Parameters
    ----------
    x, y : meshgrid arrays
        vectors or images of coordinates; should be matching numpy ndarrays

    par : list
        vector of parameters with 1 to 5 elements, defined as follows:
        par[0]: lens strength, or 'Einstein radius'
        par[1]: (optional) x-center (default = 0.0)
        par[2]: (optional) y-center (default = 0.0)
        par[3]: (optional) axis ratio (default=1.0)
        par[4]: (optional) major axis Position Angle
                in degrees c.c.w. of x axis. (default = 0.0)


    Returns
    -------
    xg, yg : gradients at the positions (x, y)


    Notes
    -----
    This routine implements an 'intermediate-axis' convention.
      Analytic forms for the SIE potential can be found in:
        Kassiola & Kovner 1993, ApJ, 417, 450
        Kormann et al. 1994, A&A, 284, 285
        Keeton & Kochanek 1998, ApJ, 495, 157
      The parameter-order convention in this routine differs from that
      of a previous IDL routine of the same name by ASB.


    Credit
    ------
    Adam S. Bolton, U of Utah, 2009

    http://www.physics.utah.edu/~bolton/python_lens_demo/

    """
    # Set parameters:
    b = np.abs(par[0]) # can't be negative!!!
    xzero = 0. if (len(par) < 2) else par[1]
    yzero = 0. if (len(par) < 3) else par[2]
    q = 1. if (len(par) < 4) else np.abs(par[3])
    phiq = 0. if (len(par) < 5) else par[4]
    eps = 0.001 # for sqrt(1/q - q) < eps, a limit expression is used.

    # Handle q > 1 gracefully:
    if (q > 1.):
        q = 1.0 / q
        phiq = phiq + 90.0

    # Go into shifted coordinats of the potential:
    phirad = np.deg2rad(phiq)
    xsie = (x-xzero) * np.cos(phirad) + (y-yzero) * np.sin(phirad)
    ysie = (y-yzero) * np.cos(phirad) - (x-xzero) * np.sin(phirad)

    # Compute potential gradient in the transformed system:
    r_ell = np.sqrt(q * xsie**2 + ysie**2 / q)
    qfact = np.sqrt(1./q - q)

    # (r_ell == 0) terms prevent divide-by-zero problems
    if qfact >= eps:
        xtg = (b/qfact) * np.arctan(qfact * xsie / (r_ell + (r_ell == 0)))
        ytg = (b/qfact) * np.arctanh(qfact * ysie / (r_ell + (r_ell == 0)))
    else:
        xtg = b * xsie / (r_ell + (r_ell == 0))
        ytg = b * ysie / (r_ell + (r_ell == 0))

    # Transform back to un-rotated system:
    xg = xtg * np.cos(phirad) - ytg * np.sin(phirad)
    yg = ytg * np.cos(phirad) + xtg * np.sin(phirad)

    # Return value:
    return xg, yg


def apply_grav_lens(image, x_cen=0, y_cen=0, r_einstein=None, eccentricity=1,
                    rotation=0):
    """
    Apply a singular isothermal ellipsoid (SIE) gravitational lens to an image

    Parameters
    ----------
    image : np.ndarray

    x_cen, y_cen : float
        [pixel] centre of the background image relative to the centre of the
        field of view

    r_einstein : float
        [pixel] Einstein radius of lens.
        If None, r_einstein = image.shape[0] // 4

    eccentricity : float
        [1..0] The ratio of semi-minor to semi-major axis for the lens

    rotation : float
        [degrees] Rotation of lens ccw from the x axis


    Returns
    -------
    lensed_image : np.ndarray


    Example
    -------
    ::

        >>> from astropy.io import fits
        >>> im = fits.getdata("my_galaxy.fits")
        >>> im2 = apply_grav_lens(im, x_cen=30, rotation=-45, eccentricity=0.5,
                                  r_einstein=300)

    """

    if r_einstein is None:
        r_einstein = image.shape[0] // 4

    shifted_image = spi.shift(image, (x_cen, y_cen))

    nx, ny = shifted_image.shape
    w = np.linspace(-nx // 2, nx // 2, nx)
    h = np.linspace(-ny // 2, ny // 2, ny)
    x, y = np.meshgrid(w,h)

    # Get the distortions from the lens
    lpar = np.asarray([r_einstein, x_cen, y_cen, eccentricity, rotation])
    xg, yg = sie_grad(x, y, lpar)

    # Pull out the pixels from the original image and place them where the lens
    #  would put them
    i = (x-xg + nx//2).astype(int)
    j = (y-yg + ny//2).astype(int)

    lensed_image = shifted_image[j.flatten(),
                                 i.flatten()].reshape(shifted_image.shape)

    return lensed_image


def elliptical(half_light_radius, plate_scale, magnitude=10, n=4,
               filter_name="Ks", normalization="total", spectrum="elliptical",
               **kwargs):
    """
    Create a extended :class:`.Source` object for a "Galaxy"

    Parameters
    ----------
    half_light_radius : float
        [arcsec]

    plate_scale : float
        [arcsec]

    magnitude : float
        [mag, mag/arcsec2]

    n : float, optional
        Power law index. Default = 4
        - n=1 for exponential (spiral),
        - n=4 for de Vaucouleurs (elliptical)

    filter_name : str, TransmissionCurve, optional
        Default is "Ks". Values can be either:
        - the name of a SimCADO filter : see optics.get_filter_set()
        - or a TransmissionCurve containing a user-defined filter

    normalization : str, optional
        ["half-light", "centre", "total"] Where the profile equals unity
        If normalization equals:
        - "half-light" : the pixels at the half-light radius have a surface
                         brightness of ``magnitude`` [mag/arcsec2]
        - "centre" : the maximum pixels have a surface brightness of
                     ``magnitude`` [mag/arcsec2]
        - "total" : the whole image has a brightness of ``magnitude`` [mag]

    spectrum : str, EmissionCurve, optional
        The spectrum to be associated with the galaxy. Values can either be:
        - the name of a SimCADO SED spectrum : see get_SED_names()
        - an EmissionCurve with a user defined spectrum


    Optional Parameters (passed to ``sersic_profile``)
    --------------------------------------------------
    ellipticity : float
        Default = 0.5

    angle : float
        [deg] Default = 30. Rotation anti-clockwise from the x-axis

    width, height : int
        [arcsec] Dimensions of the image. Default: 512*plate_scale

    x_offset, y_offset : float
        [arcsec] The distance between the centre of the profile and the centre
        of the image. Default: (dx,dy) = (0,0)


    Returns
    -------
    galaxy_src : :class:`.Source`


    See Also
    --------
    source.sersic_profile()
    optics.get_filter_set(), source.get_SED_names()
    spectral.TransmissionCurve, spectral.EmissionCurve

    """

    params = {"n"           : n,
              "ellipticity" : 0.5,
              "angle"       : 30,
              "width"       : plate_scale * 512,
              "height"      : plate_scale * 512,
              "x_offset"    : 0,
              "y_offset"    : 0}
    params.update(kwargs)

    pixular_hlr = half_light_radius / plate_scale

    im = sersic_profile(r_eff        =pixular_hlr,
                        n            =params["n"],
                        ellipticity  =params["ellipticity"],
                        angle        =params["angle"],
                        normalization=normalization,
                        width        =params["width"] /plate_scale,
                        height       =params["height"]/plate_scale,
                        x_offset     =params["x_offset"]/plate_scale,
                        y_offset     =params["y_offset"]/plate_scale)

    if isinstance(spectrum, EmissionCurve):
        lam, spec = spectrum.lam, spectrum.val
        lam, spec = scale_spectrum(lam=lam, spec=spec, mag=magnitude,
                                   filter_name=filter_name)
    elif spectrum in get_SED_names():
        lam, spec = SED(spec_type=spectrum, filter_name=filter_name,
                        magnitude=magnitude)
    else:
        print(spectrum)
        raise ValueError("Cannot understand ``spectrum``")

    galaxy_src = source_from_image(images=im, lam=lam, spectra=spec,
                                   plate_scale=plate_scale)

    return galaxy_src


def sersic_profile(r_eff=100, n=4, ellipticity=0.5, angle=30,
                   normalization="total",
                   width=1024, height=1024, x_offset=0, y_offset=0,
                   oversample=1):
    """
    Returns a 2D array with a normalised Sersic profile

    Parameters
    ----------
    r_eff : float
        [pixel] Effective (half-light) radius

    n : float
        Power law index.
        - n=1 for exponential (spiral),
        - n=4 for de Vaucouleurs (elliptical)

    ellipticity : float
        Ellipticity is defined as (a - b)/a. Default = 0.5

    angle : float
        [deg] Default = 30. Rotation anti-clockwise from the x-axis

    normalization : str, optional
        ["half-light", "centre", "total"] Where the profile equals unity
        If normalization equals:
        - "half-light" : the pixels at the half-light radius are set to 1
        - "centre" : the maximum values are set to 1
        - "total" : the image sums to 1

    width, height : int
        [pixel] Dimensions of the image

    x_offset, y_offset : float
        [pixel] The distance between the centre of the profile and the centre
        of the image

    oversample : int
        Factor of oversampling, default factor = 1. If > 1, the model is
        discretized by taking the average of an oversampled grid.


    Returns
    -------
    img : 2D array


    Notes
    -----
    Most units are in [pixel] in this function. This differs from
    :func:`.galaxy` where parameter units are in [arcsec] or [pc]

    """

    from astropy.modeling.models import Sersic2D

    # Silently cast to integer
    os_factor = np.int(oversample)

    if os_factor <= 0:
        raise ValueError("Oversampling factor must be >=1.")

    width_os = os_factor * width
    height_os = os_factor * height
    x, y = np.meshgrid(np.arange(width_os), np.arange(height_os))

    dx = 0.5 * width_os  + x_offset * os_factor
    dy = 0.5 * height_os + y_offset * os_factor

    r_eff_os = r_eff * os_factor

    mod = Sersic2D(amplitude=1, r_eff=r_eff_os, n=n, x_0=dx, y_0=dy,
                   ellip=ellipticity, theta=np.deg2rad(angle))
    img_os = mod(x, y)

    # Rebin os_factord image
    img = _rebin(img_os, os_factor)

    thresh = np.max([img[0,:].max(), img[-1,:].max(),
                     img[:,0].max(), img[:,-1].max()])
    img[img < thresh] = 0

    if "cen" in normalization.lower():
        img /= np.max(img)
    elif "tot" in normalization.lower():
        img /= np.sum(img)

    return img


def spiral_profile(r_eff, ellipticity=0.5, angle=45,
                   n_arms=2, tightness=4., arms_width=0.1,
                   central_brightness=10, normalization='total',
                   width=1024, height=1024, oversample=1,
                   **kwargs):
    """
    Creates a spiral profile with arbitary parameters

    Parameters
    ----------
     r_eff : float
        [pixel] Effective (half-light) radius

    ellipticity : float
        Ellipticity is defined as (a - b)/a. Default = 0.5

    angle : float
        [deg] Default = 45. Rotation anti-clockwise from the x-axis

    n_arms : int
        Number of spiral arms

    tightness : float
        How many times an arm crosses the major axis. Default = 4.

    arms_width : float
        An arbitary scaling factor for how think the arms should be.
        Seems to scale with central_brightness. Default = 0.1

    central_brightness : float
        An arbitary scaling factor for the strength of the central region.
        Has some connection to ars_width. Default = 10

    normalization : str, optional
        ["half-light", "centre", "total"] Where the profile equals unity
        If normalization equals:
        - "centre" : the maximum values are set to 1
        - "total" : the image sums to 1

    width, height : int, int
        [pixel] Dimensions of the image

    x_offset, y_offset : float
        [pixel] The distance between the centre of the profile and the centre
        of the image

    oversample : int
        Factor of oversampling, default factor = 1. If > 1, the model is
        discretized by taking the average of an oversampled grid.


    Optional Parameters
    -------------------
    **kwargs are passed to sersic_profile()


    Returns
    -------
    img : np.ndarray
        A 2D image of a spiral disk


    Notes
    -----
    The intensity drop-off is dictated by a sersic profile of with indes n=1,
    i.e. an exponential drop-off. This can be altered by passing the keyword
    "n=" as an optional parameter.

    Spiral structure taken from here:
    https://stackoverflow.com/questions/36095775/creating-a-spiral-structure-in-python-using-hyperbolic-tangent


    See Also
    --------
    sersic_profile()


    """

    if ellipticity >= 1.:
        raise ValueError("ellipticiy <= 1 . This is physically meaningless")

    # create a spiral
    xx, yy = np.meshgrid(np.arange(-width/2, width/2),
                         np.arange(-height/2, height/2))
    r = np.sqrt(abs(xx)**2 + abs(yy)**2)

    spiral = np.cos( n_arms * np.arctan2(xx,yy) + tightness * np.log(r**2) ) / \
             arms_width + central_brightness

    spiral[spiral < 0] = 0

    # add an exponential drop off in light intensity for the disk
    disk = sersic_profile(r_eff=r_eff, n=1, ellipticity=0, angle=0,
                          normalization=normalization, oversample=oversample,
                          width=width, height=height, **kwargs)

    img = spiral * disk
    thresh = np.max([img[0,:].max(), img[-1,:].max(),
                     img[:,0].max(), img[:,-1].max()])
    img[img < thresh] = 0

    # rotate and tilt
    ab = 1 - ellipticity
    img= spi.zoom(img, (ab, 1), order=1)
    img = spi.rotate(img, angle, order=1)

    # normalise the flux
    img[img < 0] = 0
    img = np.nan_to_num(img)
    if "cen" in normalization.lower():
        img /= np.max(img)
    elif "tot" in normalization.lower():
        img /= np.sum(img)

    return img


def spiral(half_light_radius, plate_scale, magnitude=10,
           filter_name="Ks", normalization="total", spectrum="spiral",
           **kwargs):
    """
    Create a extended :class:`.Source` object for a "Galaxy"

    Parameters
    ----------
    half_light_radius : float
        [arcsec]

    plate_scale : float
        [arcsec]

    magnitude : float
        [mag, mag/arcsec2]

    filter_name : str, TransmissionCurve, optional
        Default is "Ks". Values can be either:
        - the name of a SimCADO filter : see optics.get_filter_set()
        - or a TransmissionCurve containing a user-defined filter

    normalization : str, optional
        ["half-light", "centre", "total"] Where in the profile equals unityy
        If normalization equals:
        - "half-light" : the pixels at the half-light radius have a surface
                         brightness of ``magnitude`` [mag/arcsec2]
        - "centre" : the maximum pixels have a surface brightness of
                     ``magnitude`` [mag/arcsec2]
        - "total" : the whole image has a brightness of ``magnitude`` [mag]

    spectrum : str, EmissionCurve, optional
        The spectrum to be associated with the galaxy. Values can either be:
        - the name of a SimCADO SED spectrum : see get_SED_names()
        - an EmissionCurve with a user defined spectrum


    Optional Parameters (passed to ``spiral_profile``)
    --------------------------------------------------
    n_arms : int
        Number of spiral arms

    tightness : float
        How many times an arm crosses the major axis. Default = 4.

    arms_width : float
        An arbitary scaling factor for how think the arms should be.
        Seems to scale with central_brightness. Default = 0.1

    central_brightness : float
        An arbitary scaling factor for the strength of the central region.
        Has some connection to ars_width. Default = 10

    ellipticity : float
        Default = 0.5

    angle : float
        [deg] Default = 30. Rotation anti-clockwise from the x-axis

    n : float
         Sersic index, default = 1 (exponential disk)

    width, height : int
        [arcsec] Dimensions of the image. Default: 512*plate_scale


    Returns
    -------
    galaxy_src : scopesim.Source


    See Also
    --------
    sersic_profile(), spiral_profile()
    optics.get_filter_set(), source.get_SED_names()
    spectral.TransmissionCurve, spectral.EmissionCurve

    """

    pixular_hlr = half_light_radius / plate_scale

    params = {"n"           : 1,
              "ellipticity" : 0.5,
              "angle"       : 30,
              "width"       : pixular_hlr,
              "height"      : pixular_hlr,
              "n_arms"      : 2,
              "tightness"   : 4.,
              "arms_width"  : 0.1,
              "central_brightness" : 10}
    params.update(kwargs)


    spiral = spiral_profile(r_eff             =pixular_hlr,
                            ellipticity       =params["ellipticity"],
                            angle             =-params["angle"],
                            normalization     =normalization,
                            width             =int(2*pixular_hlr),
                            height            =int(2*pixular_hlr),
                            n_arms            =params["n_arms"],
                            tightness         =params["tightness"],
                            arms_width        =params["arms_width"],
                            central_brightness=params["central_brightness"])

    disk = sersic_profile(r_eff        =pixular_hlr,
                          n            =1,
                          ellipticity  =params["ellipticity"],
                          angle        =params["angle"],
                          normalization=normalization,
                          width        =spiral.shape[1],
                          height       =spiral.shape[0])

    thresh = np.max((disk[0,:].max(), disk[-1,:].max(),
                     disk[:,0].max(), disk[:,-1].max()))
    disk[disk < thresh] = 0

    if isinstance(spectrum, EmissionCurve):
        lam, spec = spectrum.lam, spectrum.val
        lam, spec = scale_spectrum(lam=lam, spec=spec, mag=magnitude,
                                   filter_name=filter_name)
    elif spectrum in get_SED_names():
        lam, spec = SED(spec_type=spectrum, filter_name=filter_name,
                        magnitude=magnitude)
    else:
        print(spectrum)
        raise ValueError("Cannot understand ``spectrum``")

    gal_img = (spiral + disk).T

    galaxy_src = source_from_image(images=gal_img, lam=lam, spectra=spec,
                                   plate_scale=plate_scale)

    return galaxy_src