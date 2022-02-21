import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u

import matplotlib.pyplot as plt

from scopesim.tests.mocks.py_objects.trace_list_utils import rotate, shear, \
    stretch, translate, curvify


PLOTS = False


# def rhombify(x, y, x0, y0, angle):
#     """ Attempts to make a rhombus, but this is tricky with single lines"""
#     # ..todo:: doesn't work, fix it
#     dist = ((x - x0)**2 + (y - y0)**2)**0.5
#     radius = np.min(dist)
#     # sign = np.arctan2()
#     dlen = radius * np.tan(angle / 57.29)
#     len = ((x[-1] - x[0])**2 + (y[-1] - y[0])**2)**0.5
#
#     scale = (dlen + len) / len
#     xnew = x0 + (x-x0) * scale
#     ynew = y0 + (y-y0) * scale
#
#     return xnew, ynew


def make_trace_lists(xn, yn, wmin, wmax, smin=-1.5, smax=1.5,
                     xmin=-1, xmax=1, ymin=-1, ymax=1):
    """
    Makes a list of trace line coordinates to be used in a Trace table

    Makes input for the Trace tables with the following format::

        wave    s1 ... sN   x1 ... xN   y1 ... yN
        um      arcsec      pixel       pixel

    The number of arrays returned are 3*N+1, for N trace lines. Each trace line
    must have image plane (x, y) coordinates that correspond to a single (s)
    coordinate along the slit for a given (wave) wavelength.

    Parameters
    ----------
    xn, yn : int
        Number of coordinates along each axis.
    wmin, wmax : float
        [um] Wavelength range covered by the trace
    smax, smin : float
        [arcsec] Slit coordinates that map to the trace along the minor axis
    xmin, xmax, ymin, ymax : int
        [pixel] Coordinate space to be covered by the trace

    Returns
    -------
    sl, wl, xl, yl : list of array
        - sl: a list of 1D arrays for the slit coordinate of each trace line
        - wl: [um] a single array for the wavelengths for each trace line
        - xl: [pixel] a list of 1D arrays for the image plane x coordinates
        - yl: [pixel] a list of 1D arrays for the image plane y coordinates

    """

    ws = np.linspace(wmin, wmax, max(xn, yn))   # [um]
    ss = np.linspace(smin, smax, min(xn, yn))   # [arcsec] w.r.t centre of FOV
    xs = np.linspace(xmin, xmax, xn)            # [pix] w.r.t. centre of FOV
    ys = np.linspace(ymin, ymax, yn)            # [pix] w.r.t. centre of FOV

    wl = [ws]
    sl = [np.ones(len(ws)) * s for s in ss]
    xl = [np.ones(len(ws)) * x for x in xs]
    yl = [ys] * len(xl)

    return sl, wl, xl, yl


def make_trace_table(sl, wl, xl, yl, pixel_size=0.015 * u.mm):
    """
    Generates a Trace table from the output of ``make_trace_lists``

    Expected input is:
    * wl: 1 array
    * sl, xl, yl: N arrays for each trace

    Parameters
    ----------
    sl, wl, xl, yl : list of arrays
        The output of ``make_trace_lists``
    pixel_size : quantity
        [mm] Converts pixels to mm on the image plane

    Returns
    -------

    """
    sl = np.asarray(sl)
    xl = np.asarray(xl)
    yl = np.asarray(yl)
    xn, yn = xl.shape
    wl = np.column_stack([wl] * min(xn, yn))
    if yn < xn:
        wl = wl.T

    tbl = Table(data=[wl.flatten() * u.um,
                      sl.flatten() * u.arcsec,
                      xl.flatten() * pixel_size,
                      yl.flatten() * pixel_size],
                names=["wavelength", "s", "x", "y"])

    return tbl


def trace_0(xn=3, yn=41, wmin=0.8, wmax=2.4):
    """ A basic trace straight down the middle of a 4096x4096 chip """
    sl, wl, xl, yl = make_trace_lists(xn, yn, wmin=wmin, wmax=wmax)
    for i in range(len(xl)):
        xl[i], yl[i] = stretch(xl[i], yl[i], nx=250, ny=2000)

    tbl = make_trace_table(sl, wl, xl, yl)

    return tbl


def trace_1(xn=3, yn=16, wmin=1.2, wmax=2.4,
            x0=-1250, y0=-1000, mx=-0.5, my=-0.5):
    """ A bi-directionally sheared trace in the bottom left of the chip """
    sl, wl, xl, yl = make_trace_lists(xn, yn, wmin=wmin, wmax=wmax)
    for i in range(len(xl)):
        xl[i], yl[i] = stretch(xl[i], yl[i], nx=250, ny=750)
        xl[i], yl[i] = translate(xl[i], yl[i], dx=x0, dy=y0)
        xl[i], yl[i] = shear(xl[i], yl[i], x0=x0, y0=y0, mx=mx, my=my)

    tbl = make_trace_table(sl, wl, xl, yl)

    return tbl


def trace_2(xn=3, yn=16, wmin=0.8, wmax=1.6,
            x0=-1250, y0=1000, angle=-30):
    """ A rotated trace in the top left of the chip """
    sl, wl, xl, yl = make_trace_lists(xn, yn, wmin=wmin, wmax=wmax)  # .l = list
    for i in range(len(xl)):
        xl[i], yl[i] = stretch(xl[i], yl[i], nx=250, ny=750)
        xl[i], yl[i] = translate(xl[i], yl[i], dx=x0, dy=y0)
        xl[i], yl[i] = rotate(xl[i], yl[i], x0=x0, y0=y0, angle=angle)

    tbl = make_trace_table(sl, wl, xl, yl)

    return tbl


def trace_3(xn=3, yn=16, wmin=1.4, wmax=2.0,
            x0=1250, y0=1000, x1=0, y1=1000):
    """ A bent trace in the top right of the chip """
    sl, wl, xl, yl = make_trace_lists(xn, yn, wmin=wmin, wmax=wmax)
    for i in range(len(xl)):
        xl[i], yl[i] = stretch(xl[i], yl[i], nx=250, ny=750)
        xl[i], yl[i] = translate(xl[i], yl[i], dx=x0, dy=y0)
        xl[i], yl[i] = curvify(xl[i], yl[i], x0=x1, y0=y1)

    tbl = make_trace_table(sl, wl, xl, yl)

    return tbl


def trace_4(xn=3, yn=16, wmin=0.8, wmax=1.0,
            x0=1200, y0=-300):
    """ A short straight horizontal trace in the centre right of the chip"""
    sl, wl, xl, yl = make_trace_lists(xn, yn, wmin=wmin, wmax=wmax)
    for i in range(len(xl)):
        xl[i], yl[i] = stretch(xl[i], yl[i], nx=250, ny=750)
        xl[i], yl[i] = translate(xl[i], yl[i], dx=x0, dy=y0)
        xl[i], yl[i] = rotate(xl[i], yl[i], x0=x0, y0=y0, angle=90)

    tbl = make_trace_table(sl, wl, xl, yl)

    return tbl


def trace_5(xn=3, yn=16, wmin=2.1, wmax=2.4,
            x0=1750, y0=-1750):
    """ A sheared trace extending off the chip to the bottom right"""
    sl, wl, xl, yl = make_trace_lists(xn, yn, wmin=wmin, wmax=wmax)
    for i in range(len(xl)):
        xl[i], yl[i] = stretch(xl[i], yl[i], nx=250, ny=750)
        xl[i], yl[i] = translate(xl[i], yl[i], dx=x0, dy=y0)
        xl[i], yl[i] = shear(xl[i], yl[i], x0=x0, y0=y0, mx=-0.5)

    tbl = make_trace_table(sl, wl, xl, yl)

    return tbl


def id_table(traces_ids, descriptions=None):
    """
    Makes a Ext-1 mapping table for the Trace lists FITS file

    Parameters
    ----------
    traces_ids : list
        List of indexes of Traces that should be included in the FITS file
    descriptions : list of str, optional
        Default is None. If alternative descriptions are required

    Returns
    -------
    tbl : Table

    """
    names = np.array(["description", "extension_id",
                      "aperture_id", "image_plane_id"])
    if descriptions is None:
        descriptions = ["Centre", "Sheared", "Rotated", "Curvified",
                        "90_degrees", "Off_chip"]
    descrips = [descriptions[i] for i in traces_ids]
    n = len(descrips)
    ext_ids = np.arange(2, 2+n).astype(int)
    ap_ids = np.zeros(n).astype(int)
    imp_ids = np.zeros(n).astype(int)
    data = [descrips, ext_ids, ap_ids, imp_ids]

    tbl = Table(data=data, names=names)

    return tbl


def plot_traces(traces):
    pixel_size = 0.015
    for trace in traces:
        n = (len(trace.columns) + 1) // 3
        for i, c in zip(range(n), "brg"):
            plt.scatter(trace["x"+str(i)] / pixel_size,
                        trace["y"+str(i)] / pixel_size,
                        c=trace["wavelength"], cmap="jet")

    plt.xlim(-2048, 2048)
    plt.ylim(-2048, 2048)


def plot_traces_mm(traces):
    for trace in traces:
        n = (len(trace.columns) + 1) // 3
        for i, c in zip(range(n), "brg"):
            plt.scatter(trace["x"+str(i)],
                        trace["y"+str(i)],
                        c=trace["wavelength"], cmap="jet")

    plt.xlim(-30.72, 30.72)
    plt.ylim(-30.72, 30.72)


def make_trace_hdulist(trace_ids=None, traces=None):
    """
    Makes a Trace lists FITS file

    Default Traces which can be included are:
    - 0: A basic trace straight down the middle of a 4096x4096 chip
    - 1: A bi-directionally sheared trace in the bottom left of the chip
    - 2: A rotated trace in the top left of the chip
    - 3: A bent trace in the top right of the chip
    - 4: A short straight horizontal trace in the centre right of the chip
    - 5: A sheared trace extending off the chip to the bottom right

    Parameters
    ----------
    trace_ids : list, optional
        Default is None. List of indexes of Traces that should be included in
        the FITS file
    traces : list of Tables
        Default is None. If not None, a list of Trace tables.

    Returns
    -------
    tbl : Table

    """
    if traces is None:
        traces = [trace_0(), trace_1(), trace_2(),
                  trace_3(), trace_4(), trace_5()]
    if trace_ids is None:
        trace_ids = np.arange(len(traces))
    traces = [traces[i] for i in trace_ids]

    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header["ECAT"] = 1
    pri_hdu.header["EDATA"] = 2

    id_hdu = fits.table_to_hdu(id_table(trace_ids))

    tr_hdus = [fits.table_to_hdu(trace) for trace in traces]
    for hdu in tr_hdus:
        hdu.header["WAVECOLN"] = "wavelength"
        hdu.header["SLITPOSN"] = "s"

    hdulist = fits.HDUList([pri_hdu, id_hdu] + tr_hdus)

    if PLOTS:
        plt.figure(figsize=(15, 7))
        plt.subplot(121)
        plot_traces(traces)
        plt.subplot(122)
        plot_traces_mm(traces)
        plt.show()

    return hdulist


# make_trace_hdulist()


def make_sky_header(smin=-1.5, smax=1.5):
    pass
