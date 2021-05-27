import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.modeling import models, fitting
from scipy.interpolate import InterpolatedUnivariateSpline


def rolling_median(x, n):
    """ Calculates the rolling median of a sequence for +/- n entries """
    y = [np.median(x[max(0, i-n):min(len(x), i+n+1)]) for i in range(len(x))]
    return np.array(y)


def fill_zeros(x):
    """ Fills in zeros in a sequence with the previous non-zero number """
    for i in range(1, len(x)):
        if x[i] == 0:
            x[i] = x[i-1]
    return x


def get_max_dispersion(trace_tbls, wave_min, wave_max, dwave,
                       **kwargs):
    """
    Finds the maximum distance [mm] per wavelength unit [um] along a trace

    Looks at all the trace lines (x, y) per wavelength across the slit for each
    trace, and for every trace projected onto the image plane.
    For each wavelength in the range [wave_min, wave_max] return the largest
    dx/dwave value based on all the trace projection lines.

    Parameters
    ----------
    trace_tbls : list of fits.BinTableHDU
        List of trace position tables. Units of table [um, mm, mm]
        Each table must have columns [wavelength, x0, y0, ..., xN, yN]
    wave_min, wave_max : float
        [um] minimum wavelength to look at
    dwave : float
        [um] wavelength step size

    kwargs
    ------
    x_colname, y_colname, wave_colname : str
        The name of each column for x, y, and wavelength on the image plane
        Default column names: ["x", "y", "wavelength"]
    col_number_start : int
        Default is 0. Start of the column numbering. I.e. x0, y0, s0 etc.
        If the columns start at x1, y1, etc; set ``col_number_start=1``

    Returns
    -------
    max_grad : array
        [mm/um] The maximum sensible gradient of all spectral trace projections
    waverange : array
        [um] The wavelengths corresponding to the gradients

    """

    params = {"x_colname": "x",
              "y_colname": "y",
              "wave_colname": "wavelength",
              "col_number_start": 0}
    params.update(kwargs)

    waverange = np.arange(wave_min, wave_max, dwave)
    dispersions = []
    for tbl in trace_tbls:
        if not isinstance(tbl, Table):
            tbl = Table(tbl)

        n = len([col for col in tbl.colnames if params["y_colname"] in col])
        k = params["col_number_start"]
        # .. todo: Can't use x1, etc. anymore, we have only one column x
        colnames = ["y"+str(ii) for ii in range(k, n+k)] + \
                   ["x"+str(ii) for ii in range(k, n+k)]
        for xcol in colnames:
            xpos = tbl[xcol]
            wave = tbl[params["wave_colname"]]
            if wave[0] > wave[1]:
                wave = wave[::-1]
                xpos = xpos[::-1]

            # disp is a new range [mm] derived from the trace coordinates (x, lam)
            # and the spectral resolution dwave
            mask = (waverange >= np.min(wave)) * (waverange <= np.max(wave))
            disp = np.zeros(len(waverange))
            disp[mask] = np.interp(waverange[mask], wave, xpos)
            disp /= dwave         # [mm/um] distance / wave_unit
            dispersions += [disp]

    # find the maximum dispersion for overlapping orders by using gradients
    # .. NOTE: careful of the np.abs(). Not sure if it should be here.
    grads = np.array([np.abs(np.gradient(disp)) for disp in dispersions])
    max_grad = fill_zeros(rolling_median(np.max(grads, axis=0), 15))

    # import matplotlib.pyplot as plt
    # for grad in grads:
    #     plt.plot(waverange, grad)
    # plt.scatter(waverange, max_grad)
    # plt.show()

    # max_grad is d_pos / d_wave : change in position [mm] per micron [um]
    return max_grad, waverange


def pixel_wavelength_edges(um_per_pix, waverange, wave_min, wave_max):
    """
    Get the wavelength bin edges for pixels under (a series) of spectral traces

    Returns the wavelength bin edges needed to properly project the spectrum
    according to the provided dispersion vector ``um_per_pix``

    Note: Units must be consistent, recommended [um]

    Parameters
    ----------
    um_per_pix : list, array
    waverange : list, array
    wave_min, wave_max : float

    Returns
    -------
    wave_bin_edges : array
        [um] The wavelength bin edges

    """

    wave_bin_edges = []
    wave = wave_min
    while wave < wave_max:
        wave_bin_edges += [wave]
        wave += np.interp(wave, waverange, um_per_pix)

    return np.array(wave_bin_edges)


def get_affine_parameters(coords):
    """
    Returns rotation and shear for each MTC point along a SpectralTrace

    .. note: Restrictions of this method:

       * only uses the left most coordinates for the shear
       * rotation angle is calculated using the trace extremes

    Parameters
    ----------
    coords : dict of 2D arrays
        Each dict entry ["x", "y", "s"] contains a [N, M] 2D array of
        coordinates, where:

        * N is the number of points along the slit (e.g. ~5), and
        * M is the number of positions along the trace (e.g. >100)

    Returns
    -------
    rotations : array
        [deg] Rotation angles for M positions along the Trace
    shears : array
        [deg] Shear angles for M positions along the Trace

    """

    rad2deg = 180 / np.pi
    dxs = coords["x"][-1, :] - coords["x"][0, :]
    dys = coords["y"][-1, :] - coords["y"][0, :]
    rotations = np.arctan2(dys, dxs) * rad2deg

    dxs = np.diff(coords["x"], axis=1)
    dys = np.diff(coords["y"], axis=1)
    shears = np.array([np.arctan2(dys[i], dxs[i]) for i in range(dxs.shape[0])])
    shears = np.array(list(shears.T) + [shears.T[-1]]).T
    shears = (np.average(shears, axis=0) * rad2deg) - (90 + rotations)

    return rotations, shears


def sanitize_table(tbl, invalid_value, wave_colname, x_colname, y_colname,
                   spline_order=4, ext_id=None):

    y_colnames = [col for col in tbl.colnames if y_colname in col]
    x_colnames = [col.replace(y_colname, x_colname) for col in y_colnames]

    for x_col, y_col in zip(x_colnames, y_colnames):
        wave = tbl[wave_colname].data
        x = tbl[x_col].data
        y = tbl[y_col].data

        valid = (x != invalid_value) * (y != invalid_value)
        invalid = np.invert(valid)
        if sum(invalid) == 0:
            continue

        if sum(valid) == 0:
            warnings.warn("--- Extension {} ---"
                          "All points in {} or {} were invalid. \n"
                          "THESE COLUMNS HAVE BEEN REMOVED FROM THE TABLE \n"
                          "invalid_value = {} \n"
                          "wave = {} \nx = {} \ny = {}"
                          "".format(ext_id, x_col, y_col, invalid_value,
                                    wave, x, y))
            tbl.remove_columns([x_col, y_col])
            continue

        k = spline_order
        if wave[-1] > wave[0]:
            xnew = InterpolatedUnivariateSpline(wave[valid], x[valid], k=k)
            ynew = InterpolatedUnivariateSpline(wave[valid], y[valid], k=k)
        else:
            xnew = InterpolatedUnivariateSpline(wave[valid][::-1],
                                                x[valid][::-1], k=k)
            ynew = InterpolatedUnivariateSpline(wave[valid][::-1],
                                                y[valid][::-1], k=k)

        tbl[x_col][invalid] = xnew(wave[invalid])
        tbl[y_col][invalid] = ynew(wave[invalid])

    return tbl


# ..todo: should the next three functions be combined and return a dictionary of fits?
def xilam2xy_fit(layout, params):
    '''Determine polynomial fits of FPA position

    Fits are of degree 4 as a function of slit position and wavelength.
    '''
    xi_arr = layout[params['s_colname']]
    lam_arr = layout[params['wave_colname']]
    x_arr = layout[params['x_colname']]
    y_arr = layout[params['y_colname']]

    ## Filter the lists: remove any points with x==0
    ## ..todo: this may not be necessary after sanitising the table
    #good = x != 0
    #xi = xi[good]
    #lam = lam[good]
    #x = x[good]
    #y = y[good]

    # compute the fits
    pinit_x = models.Polynomial2D(degree=4)
    pinit_y = models.Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xilam2x = fitter(pinit_x, xi_arr, lam_arr, x_arr)
    #xilam2x.name = 'xilam2x'
    xilam2y = fitter(pinit_y, xi_arr, lam_arr, y_arr)
    #xilam2y.name = 'xilam2y'
    return xilam2x, xilam2y

def xy2xilam_fit(layout, params):
    '''Determine polynomial fits of wavelength/slit position

    Fits are of degree 4 as a function of focal plane position'''

    xi_arr = layout[params['s_colname']]
    lam_arr = layout[params['wave_colname']]
    x_arr = layout[params['x_colname']]
    y_arr = layout[params['y_colname']]

    # # Filter the lists: remove any points with x==0
    # # ..todo: this may not be necessary after sanitising the table
    # good = x != 0
    # xi = xi[good]
    # lam = lam[good]
    # x = x[good]
    # y = y[good]

    pinit_xi = models.Polynomial2D(degree=4)
    pinit_lam = models.Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xy2xi = fitter(pinit_xi, x_arr, y_arr, xi_arr)
    #xy2xi.name = 'xy2xi'
    xy2lam = fitter(pinit_lam, x_arr, y_arr, lam_arr)
    #xy2lam.name = 'xy2lam'
    return xy2xi, xy2lam


def _xiy2xlam_fit(layout, params):
    '''Determine polynomial fits of wavelength/slit position

    Fits are of degree 4 as a function of focal plane position'''

    # These are helper functions to allow fitting of left/right edges
    # for the purpose of checking whether a trace is on a chip or not.

    xi_arr = layout[params['s_colname']]
    lam_arr = layout[params['wave_colname']]
    x_arr = layout[params['x_colname']]
    y_arr = layout[params['y_colname']]

    # # Filter the lists: remove any points with x==0
    # # ..todo: this may not be necessary after sanitising the table
    # good = x != 0
    # xi = xi[good]
    # lam = lam[good]
    # x = x[good]
    # y = y[good]

    pinit_x = models.Polynomial2D(degree=4)
    pinit_lam = models.Polynomial2D(degree=4)
    fitter = fitting.LinearLSQFitter()
    xiy2x = fitter(pinit_x, xi_arr, y_arr, x_arr)
    #xy2xi.name = 'xy2xi'
    xiy2lam = fitter(pinit_lam, xi_arr, y_arr, lam_arr)
    #xy2lam.name = 'xy2lam'
    return xiy2x, xiy2lam
