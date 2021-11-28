"""
OLD_simulation.py
"""

import numpy as np
from astropy.stats import sigma_clipped_stats

import scopesim.source.OLD_templates
from scopesim.source import OLD_source
from OLD_code.OLD_simcado_user_commands import UserCommands

__all__ = ["run", "snr", "check_chip_positions", "limiting_mags"]

def run(src, mode="wide", cmds=None, opt_train=None, fpa=None,
        detector_layout="small", filename=None, return_internals=False,
        filter_name=None, exptime=None, sub_pixel=False,
        **kwargs):
    """
    Run a MICADO simulation with default parameters

    Parameters
    ----------
    src : scopesim.Source
        The object of interest

    mode : str, optional
        ["wide", "zoom"] Default is "wide", for a 4mas FoV. "Zoom" -> 1.5mas

    cmds : scopesim.UserCommands, optional
        A custom set of tests_commands for the simulation. Default is None

    opt_train : scopesim.OpticalTrain, optional
        A custom optical train for the simulation. Default is None

    fpa : scopesim.Detector, optional
        A custom detector layout for the simulation. Default is None

    detector_layout : str, optional
        ["small", "centre", "full", "tiny"] Default is "small".
        "small"   - 1x 1k-detector centred in the FoV
        "tiny"    - 128 x 128 pixels centred in the FoV
        "centre"  - 1x 4k-detector centred in the FoV
        "full"    - 9x 4k-detectors

    filename : str, optional
        The filepath for where the FITS images should be saved.
        Default is None. If None, the output images are returned to the user as
        FITS format astropy.io.HDUList py_objects.

    return_internals : bool
        [False, True] Default is False. If True, the ``UserCommands``,
        ``OpticalTrain`` and ``Detector`` py_objects used in the simulation are
        returned in a tuple: ``return hdu, (cmds, opt_train, fpa)``

    filter_name : str, TransmissionCurve
        Analogous to passing INST_FILTER_TC as a keyword argument

    exptime : int, float
        [s] Analogous to passing OBS_EXPTIME as a keyword argument

    """
    pass


def check_chip_positions(filename="st.fits", x_cen=17.084, y_cen=17.084,
                         n=0.3, mode="wide"):
    """
    Creates a series of grids of stars and generates the output images

    THe number of stars in each grid corresponds to the id number of the chip
    """


    x = [-x_cen]*1 + [0]*2 + [x_cen]*3 + \
        [-x_cen]*4 + [0]*5 + [x_cen]*6 + \
        [-x_cen]*7 + [0]*8 + [x_cen]*9

    y = [-y_cen + i*n for i in range(1)] + \
        [-y_cen + i*n for i in range(2)] + \
        [-y_cen + i*n for i in range(3)] + \
        [0 + i*n for i in range(4)] + \
        [0 + i*n for i in range(5)] + \
        [0 + i*n for i in range(6)] + \
        [y_cen + i*n for i in range(7)] + \
        [y_cen + i*n for i in range(8)] + \
        [y_cen + i*n for i in range(9)]

    lam, spec = scopesim.source.OLD_templates.SED("A0V", "Ks", 15)
    src = OLD_source.Source(lam=lam, spectra=spec, x=x, y=y, ref=[0] * len(x))

    run(src, detector_layout="full", filename=filename, mode=mode)


def _make_snr_grid_fpas(filter_names=None, mmin=22, mmax=32,
                        cmds=None, **kwargs):
    """
    Makes a series of :class:`.Detector` py_objects containing a grid of stars


    Parameters
    ----------
    filter_names : list
        Which filters to use for the images. See ``scopesim.optices.get_filter_set()``

    mmin, mmax : float
        [mag] Minimum and maximum magnitudes to use for the grid of stars

    cmds : scopesim.UserCommands
        A custom set of tests_commands for building the optical train

    Optional Parameters
    -------------------
    Any Keyword-Value pairs accepted by a
    :class:`~scopesim.tests_commands.UserCommands` object

    Returns
    -------
    fpas : list
        A list of :class:`Detector` py_objects with the grid of stars for each filter
        len(fpas) == len(filter_names)
    grid : scopesim.Source
        A :class:`Source` object containing the grids of stars

    See Also
    --------
    :class:`~scopesim.tests_commands.UserCommands`

    """
    if filter_names is None:
        filter_names = ["J", "H", "Ks"]

    if isinstance(filter_names, str):
        filter_names = [filter_names]

    if not isinstance(cmds, list):
        cmds = [cmds] * len(filter_names)

    fpas = []
    grids = []
    for filt, cmd in zip(filter_names, cmds):
        if cmd is None:
            cmd = UserCommands()
        #cmd["FPA_USE_NOISE"] = "no"
        cmd["OBS_NDIT"] = 1
        cmd["FPA_LINEARITY_CURVE"] = "none"
        cmd["FPA_CHIP_LAYOUT"] = "small"
        cmd.update(kwargs)

        star_sep = cmd["SIM_PIXEL_SCALE"] * 100

        grid = scopesim.source.OLD_templates.star_grid(100, mmin, mmax, filter_name=filt, separation=star_sep)
        grids += [grid]

        hdus, (cmd, opt, fpa) = run(grid, filter_name=filt, cmds=cmd,
                                    return_internals=True)
        fpas += [fpa]

    return fpas, grid


def _get_limiting_mags(fpas, grid, exptimes, filter_names=None,
                       mmin=22, mmax=32, AB_corrs=None, limiting_sigma=5):
    """Return the limiting magnitude(s) for filter(s) and exposure time(s)


    Parameters
    ----------
    fpas : list
        The output from A list of :class:`Detector` py_objects with the grid of stars
        for each filter

    grid : scopesim.Source
        The :class:`Source` object containing the grid of stars - used for the pixel
        positions of the stars

    exptimes : array
        [s] An array of exposure times in seconds

    filter_names : list
        A list of filters. See :func:`scopesim.optics.get_filter_set`

    mmin, mmax : float
        [mag] the minimum and maximum magnitudes in the grid of stars

    AB_corrs : list
        [mag] A list of magnitude corrections to convert from Vega to AB magnitudes

    limiting_sigma : float
        [sigma] The number of sigmas to use to define the limiting magnitude.
        Default is 5*sigma


    Returns
    -------
    mags_all : list
        [mag] A list of limiting magnitudes for each exposure time for each filter
        Dimensions are [n, m] where n is the number of filters and m is the number
        of exposure times passed


    """

    from scipy import stats


    if AB_corrs is None:
        AB_corrs = np.zeros(len(fpas))
    if np.isscalar(AB_corrs):
        AB_corrs = [AB_corrs]*len(fpas)

    if filter_names is None:
        filter_names = ["J", "H", "Ks"]

    if isinstance(filter_names, str):
        filter_names = [filter_names]*len(fpas)

    if np.isscalar(exptimes):
        exptimes = [exptimes]*len(fpas)

    mags_all = []
    for fpa, filt, AB_corr in zip(fpas, filter_names, AB_corrs):

        lim_mags = []
        for exptime in exptimes:

            hdus = fpa.read_out(OBS_EXPTIME=exptime)
            im = hdus[0].data

            im_width = hdus[0].data.shape[0]
            #x = (grid.x_pix+im_width//2).astype(int)
            x = grid.x_pix.astype(int)
            y = grid.y_pix.astype(int)

            sigs, nss, snrs, bgs = [], [], [], []
            for n in range(len(x)):
                dw = 5
                w = max(dw + 5, int((1. - n/len(x)) * 20))

                ps = im[y[n]-w:y[n]+w, x[n]-w:x[n]+w]

                sig = np.copy(ps[dw:-dw, dw:-dw])
                bg = np.copy(ps)
                bg[dw:-dw, dw:-dw] = 0

                bgs += [np.average(bg[bg != 0])]
                nss += [np.std(bg[bg != 0]) * np.sqrt(np.sum(bg == 0))]
                sigs += [np.sum(sig - bgs[-1])]

            nss = np.array(nss)
            sigs = np.array(sigs)
            snr = sigs/nss

            mags = np.linspace(mmin, mmax, len(x)) + AB_corr

            mask = snr > 5
            try:
                q = stats.linregress(mags[mask], np.log10(snr[mask]))
                lim_mag = (np.log10(limiting_sigma) - q.intercept) / q.slope
            except:
                lim_mag = 0

            lim_mags += [lim_mag]
            print(exptime, filt, lim_mag)
        mags_all += [lim_mags]

    return mags_all


def plot_exptime_vs_limiting_mag(exptimes, limiting_mags,
                                 filter_names=None,
                                 colors="bgrcymk", mmin=22, mmax=29,
                                 legend_loc=3, marker="+"):
    """
    Plots exposure time versus limiting magnitudes


    Parameters
    ----------
    exptimes : list, array
        [s] Exposure times corresponding to the signal-to-noise values

    limiting_mags : array, list of array
        [mag] Limiting magnitudes for one, or more, filters for the given exposure times
        Dimensions are (1, n) for a single filter, or (m, n) for m filters

    filter_names : list
        A list of m filters. See :func:`scopesim.optics.get_filter_set`

    colors : list
        The colours to use for dots in the plot

    mmin, mmax : float
        The minimum and maximum magnitudes for the y axis

    marker : str
        The matplotlib scatter marker key

    legend_loc : int
        Location of the legend. If ``None`` is passed, no legend is plotted

    """

    import matplotlib.pyplot as plt

    # Set defaults
    if filter_names is None:
        filter_names = ["J", "H", "Ks"]

    if len(np.shape(limiting_mags)) == 1:
        limiting_mags = [limiting_mags]
    if filter_names is None:
        filter_names = ["Filter "+str(i) for i in range(np.shape(limiting_mags)[0])]

    elif isinstance(filter_names, str):
        filter_names = [filter_names]*np.shape(limiting_mags)[0]

    exptimes = np.array(exptimes)


    fig = plt.gcf()

    #ax = fig.add_axes([a_left, a_bottom, ax_width, ax_height])
    ax1 = fig.add_axes([0, 0, 1, 1])

    for mag, clr, filt in zip(limiting_mags, colors, filter_names):
        plt.plot(exptimes/3600, mag, clr+marker, label=filt)

    if legend_loc is not None:
        plt.legend(loc=legend_loc, scatterpoints=1)

    plt.xlabel("Exposure time [hours]")
    plt.ylabel("Limiting Magnitudes")
    plt.xlim(np.min(exptimes/3600) - 0.1, np.max(exptimes/3600) + 0.1)
    plt.ylim(22, 31)

    plt.grid("on")

    ax2 = fig.add_axes([0.5, 0.15, 0.45, 0.35])

    for mag, clr in zip(limiting_mags, colors):
        plt.plot(exptimes, mag, clr+marker)

    plt.plot((60 * 1, 60 * 1), (mmin, mmax), "k:")
    plt.text(60 * 1 - 5, mmin + 0.5, "1 min", horizontalalignment="right")
    plt.plot((60 * 4, 60 * 4), (mmin, mmax), "k:")
    plt.text(60 * 4 - 5, mmin + 0.5, "4 min", horizontalalignment="right")
    plt.plot((60 * 15, 60 * 15), (mmin, mmax), "k:")
    plt.text(60 * 15 - 5, mmin + 0.5, "15 min", horizontalalignment="right")

    plt.xlim(10, 1800)
    plt.ylim(mmin, mmax)
    plt.semilogx()
    plt.xlabel("Exposure time [sec]")



def limiting_mags(exptimes=None, filter_names=None,
                  AB_corrs=None, limiting_sigma=5,
                  return_mags=True, make_graph=False,
                  mmin=22, mmax=31,
                  cmds=None, **kwargs):
    """
    Return or plot a graph of the limiting magnitudes for MICADO


    Parameters
    ----------
    exptimes : array
        [s] Exposure times for which limiting magnitudes should be found

    filter_names : list
        A list of filters. See :func:`scopesim.optics.get_filter_set`

    AB_corrs : list
        [mag] A list of magnitude corrections to convert from Vega to AB magnitudes

    limiting_sigma : float
        [sigma] The number of sigmas to use to define the limiting magnitude.
        Default is 5*sigma

    return_mags : bool
        If True (defualt), the limiting magnitude are returned

    make_graph : bool
        If True (defualt), a graph of the limiting magnitudes vs exposure time is plotted
        Calls :func:`plot_exptime_vs_limiting_mag`

    cmds : scopesim.UserCommands
        A custom set of tests_commands for building the optical train


    Optional Parameters
    -------------------
    Any Keyword-Value pairs accepted by a :class:`~scopesim.UserCommands` object


    Returns
    -------
    mags_all : list
        [mag] If ``return_mags=True``, returns a list of limiting magnitudes for
        each exposure time for each filter
        Dimensions are [n, m] where n is the number of filters and m is the number
        of exposure times passed


    Notes
    -----
    Vega to AB = {"J" : 0.91 , "H" : 1.39 , "Ks" : 1.85}


    Examples
    --------
    :
        >>> # Set 30 logarithmic time bins between 1 sec and 5 hours
        >>> exptimes = np.logspace(0, np.log10(18000), num=30, endpoint=True)
        >>> limiting_mags(exptimes=exptimes, filter_names=["J", "PaBeta"],
        ...               make_graph=False)

    """
    # Set default values
    if exptimes is None:
        exptimes = [1, 60, 3600, 18000]
    if filter_names is None:
        filter_names = ["J", "H", "Ks"]

    fpas, grid = _make_snr_grid_fpas(filter_names, cmds=cmds,
                                     mmin=mmin, mmax=mmax, **kwargs)
    limiting_mags = _get_limiting_mags(fpas, grid, exptimes, filter_names,
                                       mmin=mmin, mmax=mmax, AB_corrs=AB_corrs,
                                       limiting_sigma=limiting_sigma)

    if make_graph:
        plot_exptime_vs_limiting_mag(exptimes, limiting_mags, filter_names,
                                     mmin=mmin, mmax=mmax)

    if return_mags:
        return limiting_mags



def snr_curve(exptimes, mmin=20, mmax=30, filter_name="Ks",
              aperture_radius=4, cmds=None, **kwargs):
    """
    Get the signal to noise ratios for a series of magnitudes and exposure times

    This function "observes" a grid of 100 stars equally spaced in the range [``mmin``, ``mmax``]
    The stars are observed for all times given in ``exptime`` and the SNR for each star is returned
    for each exposure time.

    Parameters
    ----------
    exptimes : float, list
        [s] exposure times for the stars

    mmin, mmax
        [mag] minimum and maximum magnitudes tor the SNR curve

    filter_name : str
        The name of a filter installed in ScopeSim - see ``:func:~scopesim.optics.get_filter_set()``

    aperture_radius : int
        [pixels] The radius of the aperture places around each star

    cmds : UserCommands object
        Used to control the observation fo the grid of stars


    Optional Parameters
    -------------------
    **kwargs : Anything that you want to pass to a ``UserCommands`` object


    Returns
    -------
    snr_array : list of arrays
        The best fit to the Magnitude-SNR curve for each entry in ``exptimes``

    mags : np.ndarray
        [mag] The magnitudes that the SNR values correspond to. Generally 100 values
        in a np.linspace range between ``mmin`` and ``mmax``


    See Also
    --------
    ``:func:~scopesim.optics.get_filter_set()``

    """

    if isinstance(exptimes, (float, int)):
        exptimes = [exptimes]

    paranal_bg = {"J" : 16.5, "H" : 14.4, "Ks" : 13.6}

    default_cmds = UserCommands()
    default_cmds["ATMO_EC"] = "none"
    default_cmds["FPA_USE_NOISE"] = "no"
    if filter_name in paranal_bg.keys():
        default_cmds["ATMO_BG_MAGNITUDE"] = paranal_bg[filter_name]

    if cmds is not None:
        default_cmds.update(cmds)

    default_cmds.update(kwargs)

    q = _make_snr_grid_fpas(filter_names=[filter_name],
                            mmin=mmin, mmax=mmax, cmds=default_cmds)
    fpa, src = q[0][0], q[1]

    mags = np.linspace(mmin, mmax, 100)

    r = aperture_radius
    r_out = 48
    r_width = 5

    snr_array = []
    for exptime in exptimes:

        hdu = fpa.read_out(OBS_EXPTIME=exptime)
        data = hdu[0].data

        sq_aps = []
        bg_stats = []

        for i in range(len(src.x_pix)):

            x, y = int(src.x_pix[i]), int(src.y_pix[i])
            sq_ap = np.copy(data[y-r:y+r+1, x-r:x+r+1])
            sq_aps += [sq_ap]

            bg_ap = np.copy(data[y-r_out:y+r_out+1, x-r_out:x+r_out+1])
            bg_ap[r_width:-r_width, r_width:-r_width] = 0

            av, med, std = sigma_clipped_stats(bg_ap[bg_ap != 0])
            bg_stats += [[av, med, std]]

        RON = default_cmds["FPA_READOUT_MEDIAN"]
        bg_med = np.array([s[1] for s in bg_stats])
        bg_std = np.array([s[0] for s in bg_stats])
        n_pix = (r*2+1)**2

        raw = np.array([np.sum(s) for s in sq_aps])
        sig = raw - bg_med * n_pix

        sig_shot = np.sqrt(sig)
        bg_shot = np.sqrt(bg_med * n_pix)
        bg_shot_std = np.sqrt(n_pix) * bg_std
        e_shot = np.sqrt(n_pix  * RON**2)

        tot_err = np.sqrt(sig_shot**2 + bg_shot**2 + e_shot**2)

        snr_val = sig / tot_err
        mask = snr > 10

        log_snr = np.log10(snr_val[mask])
        p = np.polyfit(mags[mask], log_snr, 2)
        snr_fit = 10**np.polyval(p, mags)

        snr_array += [snr_fit]

    return snr_array, mags



def plot_snr_curve(snr_array, mags, snr_markers=None):
    """
    Plots a single ``snr_curve()`` result

    Parameters
    ----------
    snr_curve : np.ndarray
        A single array from the result of ``snr_curve()`` which holds the SNR for
        each magnitude given in ``mags`` for the corresponding exposure time

    mags : np.ndarray
        [mag] An array of magnitudes generated by ``snr_curve()``

    snr_markers : list
        The SNR that should be emphasised on the graph. Default is [5,10,250]


    Example
    -------

        >>> from scopesim.simulation import snr_curve, plot_snr_curve
        >>> snrs, mag = snr_curve(exptimes=[60, 600, 3600], filter_name="J")
        >>> plot_snr_curve(snrs[-1], mag, snr_Markers=[5,10,50,250])


    """
    from matplotlib import pyplot as plt

    # Set defaults
    if snr_markers is None:
        snr_markers = [5, 10, 250]

    plt.plot(mags, snr_array, "b")
    #plt.plot(x, yfit, "k")
    plt.semilogy()

    mags_markers = np.interp(snr_markers[::-1], snr_array[::-1],
                             mags[::-1])[::-1]

    for snr_val, m, c in zip(snr_markers, mags_markers, "ryg"):
        plt.axhline(snr_val, c=c)
        plt.axvline(m, c=c)
        plt.text(mags[2], 1.1*snr_val, str(int(snr))+r"$\sigma$", color=c)
        plt.text(m - 0.1, 1.3, str(m)[:4], color=c, rotation=90,
                 verticalalignment="bottom", horizontalalignment="right")

    plt.ylim(ymin=1)

    plt.xlabel("Magnitude")
    plt.ylabel("Signal to Noise Ratio")


def plot_snr_rainbow(exptimes, mags, snr_array, snr_levels=None,
                     text_height=None):
    """
    Plot a nice rainbow curve of the SNR as a function of exposure time and magnitude

    Basically accepts the output from ''func::~scopesim.simulation.snr_curve()''

    Parameters
    ----------
    exptimes : list, np.ndarray
        Exposure times, in whatever unit you want them to be displayed

    mags : list, np.array
        [mag] A list of the magnitudes for which the SNR has been calculated

    snr_array : 2D np.ndarray
        A 2D (n,m) array where n is the length of ''exptimes'' and m is the length of ''mags''

    snr_levels : list, np.ndarray
        Which contours should be plotted on the graph. Default is [5,10,250] sigma.

    text_height : list, np.ndarray
       [mag] The height at which the contour labels should be plotted. Default is ''None''.


    Returns
    -------
    fig : matplotlib.Figure object


    """
    from matplotlib import pyplot as plt
    from matplotlib.colors import LogNorm

    # Set defaults
    if snr_levels is None:
        snr_levels = [5, 10, 250]

    fig = plt.figure(figsize=(10, 5))
    plt.contour(exptimes, mags, np.array(snr_array).T, snr_levels,
                colors=list("krygbkkkkkkkk"))

    lvls = (list(range(1, 10)) + list(range(10, 100, 10))
            + list(range(100, 1001, 100)))
    plt.contourf(exptimes, mags, np.array(snr_array).T, levels=lvls,
                 norm=LogNorm(), alpha=0.5, cmap="rainbow_r")
    clb = plt.colorbar()
    clb.set_label(r"Signal to Noise Ratio ($\sigma$)")

    if text_height is not None:
        for m, s in zip(text_height, snr_levels[::-1]):
            plt.text(3, m, str(s)+r"$\sigma$", rotation=10)

    plt.grid()
    plt.semilogx()


def mags_from_snr_array(snr_val, snr_array, mags):
    """
    Returns magnitudes that will have a certain snr from an array returned by ``snr_curve()``

    Note - this is all good, as long as you know the exposure times that
    correspond to the SNR values

    Parameters
    ----------
    snr_val : float
        The desired SNR contour

    snr_array, mags : np.ndarray
        The outputs from ``snr_curve()``

    """

    mags_fit = [np.interp(snr_val, s[::-1], mags[::-1]) for s in snr_array]

    return mags_fit


def snr(exptimes, mags, filter_name="Ks", cmds=None, **kwargs):
    """
    Returns the signal-to-noise ratio(s) for given exposure times and magnitudes

    Each time this runs, scopesim runs a full simulation on a grid of stars. Therefore
    if you are interested in the SNR for many difference expoure times and a range of
    magnitudes, it is faster to pass all of them at once to this function. See the
    example section below.

    Parameters
    ----------
    exptimes : float, list
        [s] A single or multiple exposure times

    mags : float, list
        [mag] A single or multiple magnitudes

    filter_name : str, optional
        The name of the filter to be used - See :func:`~scopesim.optics.get_filter_set`
        The default is "Ks"

    cmds : UserCommands object, optional
        Extra tests_commands to be passed to :func:`scopesim.simulation.run`.

    Optional Parameters
    -------------------
    aperture_radius
        [pixels] Default is 4. See :func:`.snr_curve`

    **kwargs : Any keyword-value pairs to be passed to the internal :class:`.UserCommands` object


    Returns
    -------
    snr_return : list
        A list of SNR values for each exposure time and each magnitude

    Examples
    --------

    A basic example of wanting the SNR for a Ks=24 star in a 1 hr observation

        >>> snr(exptimes=3600, mags=24)
        [72.69760133863036]

    However this is slow because it runs a full simulation. Hence it is better to do more at once
    If we want the SNR for the range of magnitudes J=[15, 20, 25, 30] for a 1 hr observation:

        >>> snr(exptimes=3600, mags=[15,20,25,30], filter_name="J")
        [array([  2.35125027e+04,   2.74921916e+03,   8.97552604e+01,
          8.18183097e-01])]

    Now if we were interested in different exposure times, say 10 minutes and 5 hours, for a
    24th magnitude star in the narrow band Br$\\gamma$ filter:

        >>> # Chekc the name of the Brackett Gamma filter
        >>> [name for name in scopesim.optics.get_filter_set() if "Br" in name]
        ['BrGamma']
        >>> snr(exptimes=[600, 18000], mags=24, filter_name="BrGamma")
        [8.016218764390803, 42.71569256185457]

    """


    if isinstance(exptimes, (int, float)):
        exptimes = np.array([exptimes])

    if isinstance(mags, (int, float)):
        mmin, mmax = mags - 2, mags + 2
    elif isinstance(mags, (tuple, list, np.ndarray)):
        mmin, mmax = np.min(mags), np.max(mags)
    else:
        raise ValueError("Couldn't use type(mags): "+str(type(mags)))

    snr_array, mags_array = snr_curve(exptimes, mmin=mmin, mmax=mmax,
                                      filter_name=filter_name, cmds=cmds, **kwargs)

    snr_return = []
    for i in range(len(exptimes)):
        snr_i = snr_array[i]
        snr_fit = np.interp(mags, mags_array, snr_i)
        snr_return += [snr_fit]

    return snr_return





# def snr_old(mags, filter_name="Ks", total_exptime=18000, ndit=1, cmds=None):
    # """
    # Return the signal-to-noise for a list of magnitudes in a specific filter

    # Uses the standard setup for MICADO and calculates the signal-to-noise
    # ratio or a list of magnitudes in ``mags`` in a certain broadband
    # ``filter_name``.
    # A custom UserCommands object can also be used. Note that this runs a basic
    # ScopeSim simulation len(mags) times, so execution time can be many minutes.

    # Parameters
    # ----------
    # mags : array-like
        # [vega mags] The magnitude(s) of the source(s)

    # filter_name : str, optional
        # Default is "Ks". Acceptable broadband filters are UBVRIzYJHKKs

    # exptime : float
        # [s] Total exposure time length. Default is 18000s (5 hours)

    # ndit : int, optional
        # Number of readouts during the period ``exptime``. Default is 1

    # cmds : scopesim.UserCommands, optional
        # A custom set of tests_commands for the simulations. If not specified, ScopeSim
        # uses the default MICADO parameters

    # Returns
    # -------
    # sn : np.ndarray
        # An array of signal-to-noise ratios for the magnitudes given

    # """
    # TODO: What about argument cmds? (OC)

    # logging.warning("""This is in the process of being depreciated.
                     # Use 'snr_curve()' until ScopeSim v0.5 is released""")

    # if cmds is None:
        # cmd = sim.UserCommands()
    # else:
        # cmd = cmds
    # cmd["OBS_EXPTIME"] = total_exptime / ndit
    # cmd["OBS_NDIT"] = ndit
    # cmd["INST_FILTER_TC"] = filter_name

    # opt = sim.OpticalTrain(cmd)

    # if type(mags) not in (list, tuple, np.ndarray):
        # mags = [mags]

    # sn = []
    # for mag in mags:
        # st = sim.source.star(mag)

        # fpa = sim.Detector(cmd)
        # st.apply_optical_train(opt, fpa, chips=0)
        # hdu = fpa.read_out()

        # im = hdu[0].data
        # cx, cy = np.array(im.shape) // 2
        # n = 5
        # sig = np.sum(im[cx-n:cx+n+1, cy-n:cy+n+1])
        # av = np.average(im[:200, :50])
        # std = np.std(im[:200, :50])    ## unused (OC)

        # n_pix = (2*n+1)**2
        # only_sig = sig - av*n_pix
        # only_noise = av# * np.sqrt(n_pix)  ## TODO: incorrect (OC)

        # sn += [only_sig/only_noise]

    # return np.array(sn)
