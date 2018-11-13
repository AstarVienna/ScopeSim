import warnings
import numpy as np
import pysynphot as psp
from astropy.io import ascii


def box(centre, width, waveset):
    """
    Create a Box filter curve

    Parameters
    ----------
    centre : float
        [um] centre of the box profile

    width : float
        [um] width of the box profile

    waveset : tuple
        [um] (blue, red) edges of the wavelength range. Must be outside the
        edges of the box profile

    Returns
    -------
    curve : :class:`pysynphot.ArrayBandpass`

    """

    blue_edge, red_edge = centre-width/2, centre+width/2
    if blue_edge < waveset[0]:
        warnings.warn("Blue edges outside waveset. Expanding waveset")
        waveset[0] = 0.998 * blue_edge
    if red_edge > waveset[1]:
        warnings.warn("Red edges outside waveset. Expanding waveset")
        waveset[1] = 1.002 * blue_edge

    val = np.array([0,0,1,1,0,0])
    lam = np.array([waveset[0], 0.999*blue_edge, blue_edge,
                    red_edge,   1.001*red_edge,  waveset[1]])

    curve = psp.ArrayBandpass(wave=lam, throughput=val, waveunits="um")

    return curve


def export_filter_files(fname, waveset=None, **kwargs):
    """
    Export a series of box filter curves from a list of all filters.

    Parameters
    ----------
    fname : str
        Filename of the description of the filters. Must contain the columns:
        - ``name`` - name of the filter
        - ``lam_cen`` - central wavelength in [um]
        - ``del_lam`` - width of the Box filter

    waveset : tuple
        [um] Contains the wavelength range that the saved files should contain
        in the following format (lam_min, lam_max, num_bins)

    kwargs
    ------
    Anything accepted by ``numpy.savetxt``


    Example
    -------
    ::

        >>> export_filter_files("filts_bg2.txt", waveset=(0.8,2.5))

    """

    fi = ascii.read(fname)
    filters = {fi["name"][i] : box(fi["lam_cen"][i], fi["del_lam"][i],
                                   waveset=waveset) \
               for i in range(len(fi))}

    for filt_name in filters:
        filt = filters[filt_name]
        filt_arr = np.array((filt.wave,
                             filt.throughput)).T
        print(kwargs)

        np.savetxt("TC_filter_"+filt_name+".txt", filt_arr, **kwargs)
