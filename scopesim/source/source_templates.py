from pathlib import Path

import numpy as np

from astropy import units as u
from astropy.table import Table
from astropy.io import fits
from astropy.utils.decorators import deprecated_renamed_argument

from synphot import SourceSpectrum, ConstFlux1D, Empirical1D
from synphot.units import PHOTLAM

from ..optics import image_plane_utils as ipu
from .source_utils import make_img_wcs_header
from .source import Source
from .. import rc


def empty_sky(flux=0):
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
    sky = Source(lam=[0.3, 3.0]*u.um, spectra=[flux, flux]*PHOTLAM,
                 x=[0], y=[0], ref=[0], weight=[1])
    return sky


@deprecated_renamed_argument("mag", "flux", "0.1.5")
def star(x=0, y=0, flux=0):
    """
    Source object for a single star in either vega, AB magnitudes, or Jansky

    The star is associated with the reference spectrum for each photometric
    system, therefore a reference wavelength or filter does not need to be given

    Parameters
    ----------
    x, y : float
        [arcsec] position from centre of field of view
    flux : float
        [vega mag, AB mag, Jy] Stellar brightness

    Returns
    -------
    src : Source
        A source object with a single entry table field and a reference spectrum

    """

    mag_unit = u.mag
    spec_template = vega_spectrum
    if isinstance(flux, u.Quantity):
        if flux.unit.physical_type == "spectral flux density":  # ABmag and Jy
            mag_unit = u.ABmag
            spec_template = ab_spectrum
            flux = flux.to(u.ABmag)
        flux = flux.value

    spec = spec_template()

    w = 10**(-0.4 * flux)
    ref = 0

    tbl = Table(data=[[x], [y], [w], [ref], [flux]],
                dtype=[float, float, float, int, float],
                names=["x", "y", "weight", "ref", "mag"],
                units=[u.arcsec, u.arcsec, None, None, mag_unit])
    tbl.meta["photometric_system"] = "vega" if mag_unit == u.mag else "ab"
    src = Source(spectra=spec, table=tbl)
    src.meta.update({"function_call": f"star({x=}, {y=}, {flux=})",
                     "module": "scopesim.source.source_templates",
                     "object": "star"})

    return src


def star_field(n, mmin, mmax, width, height=None, use_grid=False):
    """
    Creates a super basic field of stars with random positions and brightnesses

    Parameters
    ----------
    n : int
        number of stars

    mmin, mmax : float, astropy.Quantity
        [mag, ABmag, Jy] min and max magnitudes/fluxes of the population stars.
        If floats, then assumed Quantity is vega magnitudes

    width : float
        [arcsec] width of region to put stars in

    height : float, optional
        [arcsec] if None, then height=width

    use_grid : bool, optional
        Place stars randomly or on a grid

    Returns
    -------
    stars : scopesim.Source object
        A Source object with a field of stars that can be fed into the method:
        OpticalTrain.observe()

    See Also
    --------
    OpticalTrain.observe
    OpticalTrain.readout

    """
    if height is None:
        height = width

    mag_unit = u.mag
    spec_template = vega_spectrum
    if isinstance(mmin, u.Quantity):
        if mmin.unit.physical_type == "spectral flux density":  # ABmag and Jy
            mag_unit = u.ABmag
            spec_template = ab_spectrum
            mmin, mmax = mmin.to(u.ABmag), mmax.to(u.ABmag)
        mmin, mmax = mmin.value, mmax.value

    spec = spec_template()

    if rc.__config__["!SIM.random.seed"] is not None:
        np.random.seed(rc.__config__["!SIM.random.seed"])

    if use_grid:
        nw = np.ceil(n**0.5)
        nh = np.ceil(n / nw)
        x, y = np.mgrid[0:1:1/nw, 0:1:1/nh] - 0.5
        positions = x.flatten(), y.flatten()
    else:
        positions = np.random.random(size=(2, n)) - 0.5

    x = width * positions[0]
    y = height * positions[1]
    mags = np.linspace(mmin, mmax, n)
    w = 10**(-0.4 * mags)
    ref = np.zeros(n, dtype=int)

    tbl = Table(data=[x, y, w, ref, mags],
                names=["x", "y", "weight", "ref", "mag"],
                units=[u.arcsec, u.arcsec, None, None, mag_unit])
    tbl.meta["photometric_system"] = "vega" if mag_unit == u.mag else "ab"
    stars = Source(spectra=spec, table=tbl)

    return stars


def uniform_illumination(xs, ys, pixel_scale, flux=None, spectrum=None):
    """
    Return a Source for a uniformly illuminated area

    Parameters
    ----------
    xs, ys : list of float
        [arcsec] min and max extent of each dimension relative to FOV centre
        E.g. `xs=[-1, 1], ys=[5, 5.5]`

    pixel_scale : float
        [arcsec]

    flux : astropy.Quantity
        [mag, ABMag, Jy] Flux per arcsecond of the Source

    Returns
    -------
    src : scopesim.Source

    Examples
    --------
    A 200x200 uniform illumination Source at 1 Jy/arcsec2
    ::

        src = uniform_illumination(xs=[-1,1], ys=[-1, 1],
                                   pixel_scale=0.01, flux=1*u.Jy)


    A source that extends just past the MICADO 15" slit dimensions with a flux
    of 10 mag/arcsec2
    ::

        src = uniform_illumination(xs=[-8, 8], ys=[-0.03, 0.03],
                                   pixel_scale=0.004, flux=10*u.mag)

    Using a self made frequency-comb spectrum with 1 Jy lines ever 0.1µm
    ::

        import numpy as np
        from astropy import units as u
        from synphot import SourceSpectrum, Empirical1D

        wave = np.arange(0.7, 2.5, 0.001) * u.um
        flux = np.zeros(len(wave))
        flux[::100] = 1 * u.Jy
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

        src = uniform_illumination(xs=[-8, 8], ys=[-0.03, 0.03],
                                   pixel_scale=0.004, spectrum=spec)

    """
    if flux is not None:
        mag_unit = u.mag
        spec_template = vega_spectrum
        if isinstance(flux, u.Quantity):
            if flux.unit.physical_type == "spectral flux density":  # ABmag and Jy
                mag_unit = u.ABmag
                spec_template = ab_spectrum
                flux = flux.to(u.ABmag)
            flux = flux.value

        spec = spec_template()
        scale_factor = pixel_scale ** 2 * 10 ** (-0.4 * flux)
    elif spectrum is not None:
        spec = spectrum
        mag_unit = "PHOTLAM"
        scale_factor = pixel_scale ** 2
    else:
        raise ValueError(f"Either flux or spectrum must be passed: {flux}, {spectrum}")

    hdr = ipu.header_from_list_of_xy(x=np.array(xs) / 3600.,
                                     y=np.array(ys) / 3600.,
                                     pixel_scale=pixel_scale / 3600.)

    data = scale_factor * np.ones([max(1, hdr["NAXIS2"]), max(1, hdr["NAXIS1"])])
    hdu = fits.ImageHDU(header=hdr, data=data)

    hdu.header["SPEC_REF"] = 0
    hdu.header["SPECMAG"] = 0
    hdu.header["SPECUNIT"] = str(mag_unit)

    src = Source(image_hdu=hdu, spectra=[spec])

    return src


def uniform_source(sp=None, extent=60):
    """
    Simplified form of scopesim_templates.misc.uniform_source, mostly intended for testing

    This function creates an image with extend^2 pixels with pixel size of 1 arcsec^2 so provided amplitudes
    are in flux or magnitudes per arcsec^2

    It accepts any synphot.SourceSpectrum compatible object

    sp : synphot.SourceSpectrum
         defaults to vega_spectrum() with magnitude 0 mag/arcsec2

    extent : int, default 60
        extension of the field in arcsec, will always produce a square field. Default value produces a field of 60x60 arcsec

    """
    if sp is None:
        sp = vega_spectrum()

    data = np.ones(shape=(extent, extent))

    header = make_img_wcs_header(pixel_scale=1, image_size=data.shape)
    hdu = fits.ImageHDU(header=header, data=data)

    src = Source(spectra=sp, image_hdu=hdu)

    return src


def vega_spectrum(mag=0):
    if isinstance(mag, u.Quantity):
        mag = mag.value
    # HACK: Turn Path object back into string, because not everything
    #       that depends on this function can handle Path objects (yet)
    vega = SourceSpectrum.from_file(str(Path(rc.__pkg_dir__, "vega.fits")))
    vega = vega * 10 ** (-0.4 * mag)
    return vega


def st_spectrum(mag=0):
    waves = np.geomspace(100, 300000, 50000)
    sp = ConstFlux1D(amplitude=mag*u.STmag)

    return SourceSpectrum(Empirical1D, points=waves, lookup_table=sp(waves))


def ab_spectrum(mag=0):
    waves = np.geomspace(100, 300000, 50000)
    sp = ConstFlux1D(amplitude=mag * u.ABmag)

    return SourceSpectrum(Empirical1D, points=waves, lookup_table=sp(waves))
