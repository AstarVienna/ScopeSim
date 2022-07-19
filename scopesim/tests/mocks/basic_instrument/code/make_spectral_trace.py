"""
3 traces vertically covering the regions:
- J: 1.0, 1.35     dwave = 0.35 --> R~2800
- H: 1.3, 1.8      dwave = 0.5  --> R~2000
- K: 1.75, 2.5     dwave = 0.75 --> R~1300

The detector is 1024 * 0.2 = 204.8" [+/-102.4"] and 10.24mm [+/-5.12mm]
Slit is 60", or 300 pixel, so trace x-borders should be:
- J: (-95, -35) ["], --> /20 "/mm --> (-4.75, -1.75) [mm] centred at -3.25
- H: (-30,  30) ["], --> /20 "/mm --> (-1.5,   1.5)  [mm] centred at 0
- K: ( 35,  95) ["], --> /20 "/mm --> ( 1.75,  4.75) [mm] centred at 3.25

Oliver's SpectralTraces follow a O(5) polynomial, therefore we need 5
points along the slits
s = [-30, -10, 0, 10, 30] ["]
dx = [-1.5, -0.5, 0, 0.5, 1.5 ] [mm]

The dispersion should follow a 10-point logspace increase from
wave_min, wave_max for each of the trace extremes

y value go *linspace* from +/-500 px --> +/-5 mm

The format of the fits file is:
wave, slit ["], x [mm], y [mm]

"""
from copy import deepcopy
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from matplotlib import pyplot as plt


def make_lss_trace_file():
    names = ["J", "H", "K"]
    wave_mins = np.array([1.0, 1.3, 1.75])
    wave_maxs = np.array([1.35, 1.8, 2.5])
    s_cens = np.array([-65., 0., 65.])
    dss = np.array([-30., 10., 0., 10., 30.])

    hdus = []
    for wave_min, wave_max, s_cen, name in zip(wave_mins, wave_maxs, s_cens, names):
        n = 101

        waves = np.geomspace(wave_min, wave_max, n)
        waves = np.array([[w]*5 for w in waves]).flatten()

        ys = np.array([[y]*5 for y in np.linspace(-5, 5, n)]).flatten()

        ss = np.array(list(dss) * n)
        xs = (ss + s_cen) / 20.

        tbl = Table(names=["wavelength", "s", "x", "y"],
                    data=[waves, ss, xs, ys],
                    units=[u.um, u.arcsec, u.mm, u.mm])
        hdu = fits.table_to_hdu(tbl)
        hdu.header["EXTNAME"] = f"TRACE_{name}"
        hdus += [hdu]

    tbl = Table(names=["description", "extension_id", "aperture_id", "image_plane_id"],
                data=[[hdu.header["EXTNAME"] for hdu in hdus],
                      [2,3,4],
                      [0,0,0],
                      [0,0,0]])
    cat_hdu = fits.table_to_hdu(tbl)
    cat_hdu.header["EXTNAME"] = "TOC"

    pri_hdr = fits.PrimaryHDU()
    pri_hdr.header["ECAT"] = 1         # where to find the catalogue table
    pri_hdr.header["EDATA"] = 2        # where the data tables start

    hdul = fits.HDUList([pri_hdr, cat_hdu] + hdus)
    hdul.writeto("../INS_lss_traces.fits", overwrite=True)

# make_lss_trace_file()

def make_ifu_trace_file():
    names = [f"Ap{i}" for i in range(5)]
    wave_min = 1.75
    wave_max = 2.5
    s_cens = np.array([-64., -32., 0., 32., 64.])   # across detector
    dss = np.array([-14., -7., 0., 7., 14.])     # across slit

    hdus = []
    for s_cen, name in zip(s_cens, names):
        n = 101

        waves = np.geomspace(wave_min, wave_max, n)
        waves = np.array([[w] * 5 for w in waves]).flatten()    # 5 pos along slit

        ys = np.array([[y] * 5 for y in np.linspace(-5, 5, n)]).flatten()

        ss = np.array(list(dss) * n)
        xs = (ss + s_cen) / 20.

        tbl = Table(names=["wavelength", "s", "x", "y"],
                    data=[waves, ss, xs, ys],
                    units=[u.um, u.arcsec, u.mm, u.mm])
        hdu = fits.table_to_hdu(tbl)
        hdu.header["EXTNAME"] = f"TRACE_{name}"
        hdus += [hdu]

    tbl = Table(
        names=["description", "extension_id", "aperture_id", "image_plane_id"],
        data=[[hdu.header["EXTNAME"] for hdu in hdus],
              [2, 3, 4, 5, 6],
              [0, 1, 2, 3, 4],
              [0, 0, 0, 0, 0]])
    cat_hdu = fits.table_to_hdu(tbl)
    cat_hdu.header["EXTNAME"] = "TOC"

    pri_hdr = fits.PrimaryHDU()
    pri_hdr.header["ECAT"] = 1  # where to find the catalogue table
    pri_hdr.header["EDATA"] = 2  # where the data tables start

    hdul = fits.HDUList([pri_hdr, cat_hdu] + hdus)
    hdul.writeto("../INS_ifu_traces.fits", overwrite=True)

# make_ifu_trace_file()
