import os
import glob

import numpy as np

from OLD_code import OLD_spectral as sc
from scopesim.utils import find_file
from scopesim import rc


def get_filter_curve(filter_name):
    """
    Return a Vis/NIR broadband filter TransmissionCurve object

    Parameters
    ----------
    filter_name : str

    Notes
    -----
    Acceptable filters can be found be calling get_filter_set()

    To access the values use TransmissionCurve.lam and TransmissionCurve.val

    Examples
    --------
        >>> transmission_curve = get_filter_curve("TC_filter_Ks.dat")
        >>> wavelength   = transmission_curve.lam
        >>> transmission = transmission_curve.val
    """

    fname = find_file(filter_name, silent=True)
    if fname is None:
        fname = find_file("TC_filter_" + filter_name + ".dat")
        if fname is None:
            raise ValueError("filter not recognised: " + filter_name)

    return sc.TransmissionCurve(filename=fname)


def get_filter_set(path=None):
    """
    Return a list of the filters installed in the package directory
    """
    if path is None:
        path = os.path.join(rc.__data_dir__, "data")
    lst = [i.replace(".dat", "").split("TC_filter_")[-1]
           for i in glob.glob(os.path.join(path, "TC_filter*.dat"))]
    return lst


def plot_filter_set(path=None, filters="All", cmap="rainbow", filename=None,
                    show=True):
    """
    Plot a filter transmision curve or transmision curve for a list of filters

    Parameters
    ----------
    path : str
        the location of the filters, set to None to use the default one, passed
        to get_filter_set
    filters : str or list
        a filter or a list of filters to be plotted, acceptable filters can be
        found calling get_filter_set()
    cmap : str
        any ``matplotlib`` colormap, defaulted to rainbow
    filename : str
        a filename to save the figure if necessary
    show : boolean
        if True, the plot is shown immediately in an interactive session

    Returns
    -------
    a matplotlib object

    Notes
    -----

    Examples
    --------

        >>> plot_filter_set()
        >>> plot_filter_set(cmap="viridis")
        >>> plot_filter_set(filters="Ks")
        >>> plot_filter_set(filters=("U","PaBeta","Ks"),savefig="filters.png")

    """
    import matplotlib.pyplot as plt
    from matplotlib import rcParams, cycler
    from matplotlib.cm import get_cmap

    filter_names = []
    if np.size(filters) == 1:
        filter_names = [filters]
    if np.size(filters) > 1:
        filter_names = filters

    if filters.lower() == "all":
        filter_names = get_filter_set(path)

    cmap = get_cmap(cmap)
    rcParams['axes.prop_cycle'] = \
        cycler(color=cmap(np.linspace(0, 1, np.size(filter_names))))

    peaks = np.zeros(np.size(filter_names))
    i = 0
    for filter_name in filter_names:

        tcurve = get_filter_curve(filter_name)
        wave = tcurve.lam[tcurve.val > 0.02]
        tran = tcurve.val[tcurve.val > 0.02]

        lam_peak = wave[tran == np.max(tran)]
        peaks[i] = lam_peak[0]
        i += 1

    ordered_names = [x for _, x in sorted(zip(peaks, filter_names))]

    for filter_name in ordered_names:
        tcurve = get_filter_curve(filter_name)
        wave = tcurve.lam[tcurve.val > 0.02]
        tran = tcurve.val[tcurve.val > 0.02]
        lmin = np.min(wave)
        lmax = np.max(wave)
        lam_peak = wave[tran == np.max(tran)]
        if (lmax-lmin)/lam_peak[0] > 0.1:
            plt.plot(wave, tran, "--", label=filter_name)
        else:
            plt.plot(wave, tran, "-", label=filter_name)

    plt.xlabel(r"wavelength [$\mu$m]")
    plt.ylabel("transmission")
    lgd = plt.legend(loc=(1.03, 0))
    if filename is not None:
        plt.savefig(filename, bbox_extra_artists=(lgd,), bbox_inches='tight')
    if show:
        plt.show()


def get_filter_table(path=None, filters="all"):

    """
    Return an astropy.table for a list of filters

    Parameters
    ----------
    path : str
        the location of the filters, set to None to use the default one, passed
        to get_filter_set
    filters : str or list
        a filter or a list of filters to be plotted, acceptable filters can be
        found calling get_filter_set()

    Returns
    -------
    an astropy.table

    Notes
    -----
    It will ONLY return values for filters that follow scopesim format

    Examples
    --------

    Obtaining table for a set of filters::

        >>> from OLD_code import OLD_optics_utils
        >>> filter_names = ["J", "Ks", "PaBeta", "U", "Br-gamma"]
        >>> table = OLD_optics_utils.get_filter_table(filters=filter_names)
        >>> filter_centers = table["center"].data
        >>> print(filter_centers)
        [1.24794438 2.14487698 2.16986118]

    Notice that only three values are printed as the U filter does not follow
    (yet) the scopesim format

    """

    # Obtaining format of the table
    filter_table = get_filter_curve("Ks").filter_table()
    filter_table.remove_row(0)

    filter_names = []
    if np.size(filters) == 1:
        filter_names = [filters]
    if np.size(filters) > 1:
        filter_names = filters
    if filters == "all":
        filter_names = get_filter_set(path)

    for name in filter_names:
        try:
            tcurve = get_filter_curve(name)
            table = tcurve.filter_table()
            filter_table.add_row(table[0])
        except ValueError:
            pass

    return filter_table


