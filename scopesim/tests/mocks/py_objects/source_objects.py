import numpy as np
from astropy import units as u, wcs
from astropy.io import fits
from astropy.table import Table
from synphot import SourceSpectrum, Empirical1D

from scopesim.source.source import Source
from scopesim.source.source_templates import vega_spectrum

from . import FILES_PATH


def _table_source():
    n = 101
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=4 * np.ones(n) * unit),
             SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit),
             SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n)[::-1] * unit)]
    tbl = Table(names=["x", "y", "ref", "weight"],
                data=[[5,  0, -5,  0]*u.arcsec,
                      [5, -9, 5,  0] * u.arcsec,
                      [2,  0,  1,  0],
                      [1,  1,  1,  2]])
    tbl_source = Source(table=tbl, spectra=specs)

    return tbl_source


def _table_source_overlapping():
    """A table with sources that are exactly at the same place.

    This allows testing whether the fluxes are correctly stacked.

    Four sources are stacked at (-1, -1), and one source with 4 times
    the weight is at (1, 1). Both regions should have the same flux.
    """
    n = 101
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.1150, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=1e9 * np.ones(n) * unit)]
    w1 = 0.0005
    tbl = Table(names=["x", "y", "ref", "weight", "spec_types"],
                data=[[-20, -20, -20, -20, 20] * u.arcsec,
                      [-20, -20, -20, -20, 20] * u.arcsec,
                      [0,  0,  0,  0, 0],
                      [w1,  w1,  w1, w1, 4 * w1],
                      ["F0II",] * 5,
                ])
    tbl_source = Source(table=tbl, spectra=specs)
    return tbl_source


def _basic_img_hdu(weight):
    n = 51
    im_wcs = wcs.WCS(naxis=2)
    im_wcs.wcs.cunit = [u.arcsec, u.arcsec]
    im_wcs.wcs.cdelt = [0.2, 0.2]
    im_wcs.wcs.crval = [0, 0]
    im_wcs.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    # im_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    im_wcs.wcs.ctype = ["LINEAR", "LINEAR"]

    im = np.random.random(size=(n, n)) * 1e-9 * weight
    im[n-1, 1] += 5 * weight
    im[1, 1] += 5 * weight
    im[n // 2, n // 2] += 10 * weight
    im[n // 2, n-1] += 5 * weight

    im_hdu = fits.ImageHDU(data=im, header=im_wcs.to_header())
    im_hdu.header["SPEC_REF"] = 0

    return im_hdu


def _image_source(dx=0, dy=0, angle=0, weight=1):
    """
    Produce a source with 3 point sources on a random BG.

    Parameters
    ----------
    dx, dy : float
        [arcsec] Offset from optical axis
    angle : float
        [deg]
    weight : float

    Returns
    -------
    source
    """
    n = 101
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit)]

    im_hdu = _basic_img_hdu(weight)
    im_source = Source(image_hdu=im_hdu, spectra=specs)

    angle = angle * np.pi / 180
    im_source.fields[0].header["CRVAL1"] += dx / 3600
    im_source.fields[0].header["CRVAL2"] += dy / 3600
    im_source.fields[0].header["PC1_1"] = np.cos(angle)
    im_source.fields[0].header["PC1_2"] = np.sin(angle)
    im_source.fields[0].header["PC2_1"] = -np.sin(angle)
    im_source.fields[0].header["PC2_2"] = np.cos(angle)

    return im_source


def _fits_image_source():
    n = 50
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit)]

    hdulist = fits.open(FILES_PATH / "test_image.fits")
    fits_src = Source(image_hdu=hdulist[0], spectra=specs)

    return fits_src


def _cube_source(**kwargs):
    """
    Produce a source with 3 point sources on a random BG.

    Parameters
    ----------
    dx, dy : float
        [arcsec] Offset from optical axis
    angle : float
        [deg]
    weight : float

    Returns
    -------
    source
    """
    n = 101
    im_hdu = _basic_img_hdu(kwargs.get("weight", 1))
    # im_src = _image_source(**kwargs)
    # data = im_src.fields[0].data
    data = im_hdu.data

    # taken from Source
    for i in [1, 2]:
        unit = u.Unit(im_hdu.header[f"CUNIT{i}"].lower())
        val = float(im_hdu.header[f"CDELT{i}"])
        im_hdu.header[f"CUNIT{i}"] = "deg"
        im_hdu.header[f"CDELT{i}"] = val * unit.to(u.deg)

    dx = kwargs.get("dx", 0)
    dy = kwargs.get("dy", 0)
    angle = kwargs.get("angle", 0)

    angle = angle * np.pi / 180
    im_hdu.header["CRVAL1"] += dx / 3600
    im_hdu.header["CRVAL2"] += dy / 3600
    im_hdu.header["PC1_1"] = np.cos(angle)
    im_hdu.header["PC1_2"] = np.sin(angle)
    im_hdu.header["PC2_1"] = -np.sin(angle)
    im_hdu.header["PC2_2"] = np.cos(angle)

    # Broadcast the array onto a 3rd dimension and scale along the new axis
    im_hdu.data = data[None, :, :] * np.linspace(0, 4, n)[:, None, None]
    # im_src.spectra = {}

    # FIXME: CRPIX might be wrong here, aka off-by-one!!
    # But all other code assumes it like this, so I'm keeping it for now.
    # astropy WCS spectral would need 51 to work correctly...
    cube_hdr_dict = {"CUNIT3": "um", "CTYPE3": "WAVE", "CDELT3": 0.02,
                     "CRVAL3": 1.5, "CRPIX3": 50, "SPEC_REF": None,
                     "BUNIT": "ph s-1 m-2 um-1"}

    im_hdu.header.update(cube_hdr_dict)
    im_src = Source(cube=im_hdu)

    return im_src


def _combined_source(im_angle=0, dx=(0, 0, 0), dy=(0, 0, 0), weight=(1, 1, 1)):
    tblsrc1 = _table_source()

    tblsrc2 = _table_source()
    tblsrc2.fields[0]["x"] += dx[0]
    tblsrc2.fields[0]["y"] += dy[0]
    tblsrc2.fields[0]["weight"] *= weight[0]

    tblsrc3 = _table_source()
    tblsrc3.fields[0]["x"] += dx[1]
    tblsrc3.fields[0]["y"] += dy[1]
    tblsrc3.fields[0]["weight"] *= weight[1]

    imsrc = _image_source(dx[2], dy[2], im_angle, weight[2])

    src = tblsrc1 + tblsrc2 + tblsrc3 + imsrc

    return src


def _single_table_source(weight=1, n=3):
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.ones(n) * weight * unit)]
    tbl = Table(names=["x", "y", "ref", "weight"],
                data=[[0]*u.arcsec, [0]*u.arcsec, [0], [1]])
    tbl_source = Source(table=tbl, spectra=specs)

    return tbl_source


def _unity_source(dx=0, dy=0, angle=0, weight=1, n=100):
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.ones(n) * unit)]

    im_wcs = wcs.WCS(naxis=2)
    im_wcs.wcs.cunit = [u.arcsec, u.arcsec]
    im_wcs.wcs.cdelt = [1, 1]
    im_wcs.wcs.crval = [0, 0]
    im_wcs.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    # im_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    im_wcs.wcs.ctype = ["LINEAR", "LINEAR"]

    im = np.ones((n, n))

    im_hdu = fits.ImageHDU(data=im, header=im_wcs.to_header())
    im_hdu.header["SPEC_REF"] = 0
    im_source = Source(image_hdu=im_hdu, spectra=specs)

    angle = angle * np.pi / 180
    im_source.fields[0].header["CRVAL1"] += dx * u.arcsec.to(u.deg)
    im_source.fields[0].header["CRVAL2"] += dy * u.arcsec.to(u.deg)
    im_source.fields[0].header["PC1_1"] = np.cos(angle)
    im_source.fields[0].header["PC1_2"] = np.sin(angle)
    im_source.fields[0].header["PC2_1"] = -np.sin(angle)
    im_source.fields[0].header["PC2_2"] = np.cos(angle)
    im_source.fields[0].data *= weight

    return im_source


def _empty_sky():
    n = 3
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.zeros(n) * unit)]
    tbl = Table(names=["x", "y", "ref", "weight"],
                data=[[0], [0], [0], [0]])
    tbl_source = Source(table=tbl, spectra=specs)

    return tbl_source


def _vega_source(mag=0, x=0, y=0):
    specs = [vega_spectrum(mag)]
    tbl = Table(names=["x", "y", "ref", "weight"],
                data=[[x]*u.arcsec, [y]*u.arcsec, [0], [1]])
    tbl_source = Source(table=tbl, spectra=specs)

    return tbl_source
