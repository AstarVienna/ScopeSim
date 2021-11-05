import numpy as np
from astropy.io import fits
import astropy.units as u
from photutils import aperture_photometry, CircularAperture

from synphot import SourceSpectrum, Empirical1D
import scopesim as sim
from scopesim import rc



def point_source_photometry_versus_cdelt_miguel(spline_order, instrument="METIS",
                                                instrument_mode="img_lm"):
    """
    Check dependence of point source photometry on Source pixel scale

    The function runs a sequence of simulations for varying CDELT1/2 in the ImageHDU that
    defines the Source object. The point source is defined by a single bright pixel in the input
    image whose pixel scale is given by the CDELT1/2 keywords in the associated FITS header. Ideally,
    the simulated ImagePlane should be independent of the input pixel scale. However, the spline interpolation
    used for rescaling the source image to the image plane introduces some numerical noise in the photometry of
    the simulated image. In the worst case, for linear spline interpolation and a source pixel scale that
    is less than half the image plane pixel scale, the object may vanish entirely.
    """
    # Create the instrument, explicitely set spline order to be used
    cmd = sim.UserCommands(use_instrument=instrument,
                           set_modes=[instrument_mode])
    cmd["!SIM.computing.spline_order"] = spline_order
    opttrain = sim.OpticalTrain(cmd)
    opttrain["detector_linearity"].include = False

    # Create an image of a point source - a single bright pixel. Header has arbitrary CDELT1/2 to start with.
    img = np.zeros((21, 21), dtype=np.float32)
    img[10, 10] = 1.
    header = fits.Header()
    header.update({'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CDELT1': 0.1,
                   'CDELT2': 0.1,
                   'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec', 'CRPIX1': 10,
                   'CRPIX2': 10,
                   'CRVAL1': 0., 'CRVAL2': 0.,
                   'BUNIT': 'Jy'})
    img_hdu = fits.ImageHDU(data=img, header=header)

    # Create a spectrum to assign 1 Jy to image pixel value 1
    wave = np.linspace(0.5, 18.5, 1001)
    flux = sim.source.source_templates.ab_spectrum(
        mag=(1 * u.Jy).to(u.ABmag).value)(wave * u.um)
    spec = SourceSpectrum(Empirical1D, points=wave * u.um, lookup_table=flux)
    src = sim.Source(spectra=[spec], image_hdu=img_hdu)

    # Instantiate lists to hold results
    bglevel = []
    bgnoise = []
    starsum = []
    starnoise = []
    starmax = []

    # Test for CDELT1/2 between 0.001 and 0.009 arcsec
    cdeltarr = np.arange(0.001, 0.009, 0.0002)
    aperture = CircularAperture([(1028., 1028.)], r=10.)
    # Run the loop

    cmd = sim.UserCommands(use_instrument=instrument,
                           set_modes=[instrument_mode])
    cmd["!SIM.computing.spline_order"] = spline_order
    opttrain = sim.OpticalTrain(cmd)
    opttrain["detector_linearity"].include = False
    i = 0
    for cdelt in cdeltarr:
        # Set CDELT1/2
        src.fields[0].header.update({"CDELT1": cdelt, "CDELT2": cdelt})

        t1 = time()
        opttrain.observe(src, update=True)
        t2 = time()
        # PRINTIN STUFF
        #    print("Elapsed Time:", t2 - t1)
        print("iteration:", i)
        print("opttrain._last_source:")
        print(opttrain._last_source)
        opttrain._last_source.fields[0].writeto("file_" + str(i) + ".fits")

        result = opttrain.image_planes[0].data

        starmax.append(result.max())
        bglevel.append(np.median(result[200:800, 200:800]))
        bgnoise.append(np.std(result[200:300, 200:300]))

        phot_table = aperture_photometry(result - bglevel[-1], aperture,
                                         error=np.ones_like(result) * bgnoise[
                                             -1])

        starsum.append(phot_table['aperture_sum'][0])
        starnoise.append(phot_table['aperture_sum_err'][0])
        #    print(starsum, starnoise, bglevel, bgnoise, starmax
        print("*" * 20)
        # pbar.value += 1
        i = i + 1

    return {"starsum": starsum, "starnoise": starnoise, "bglevel": bglevel,
            "bgnoise": bgnoise,
            "maxpixel": starmax}


def point_source_photometry_versus_cdelt_oliver(spline_order, instrument="METIS",
                                                instrument_mode="img_lm"):
    """
    Check dependence of point source photometry on Source pixel scale

    The function runs a sequence of simulations for varying CDELT1/2 in the ImageHDU that
    defines the Source object. The point source is defined by a single bright pixel in the input
    image whose pixel scale is given by the CDELT1/2 keywords in the associated FITS header. Ideally,
    the simulated ImagePlane should be independent of the input pixel scale. However, the spline interpolation
    used for rescaling the source image to the image plane introduces some numerical noise in the photometry of
    the simulated image. In the worst case, for linear spline interpolation and a source pixel scale that
    is less than half the image plane pixel scale, the object may vanish entirely.
    """

    # Create a spectrum to assign 1 Jy to image pixel value 1
    wave = np.linspace(0.5, 18.5, 1001)
    flux = sim.source.source_templates.ab_spectrum(
        mag=(1 * u.Jy).to(u.ABmag).value)(wave * u.um)
    spec = SourceSpectrum(Empirical1D, points=wave * u.um, lookup_table=flux)

    # Instantiate lists to hold results
    bglevel = []
    bgnoise = []
    starsum = []
    starnoise = []
    starmax = []

    # Test for CDELT1/2 between 0.001 and 0.009 arcsec
    cdeltarr = np.arange(0.001, 0.009, 0.0002)
    aperture = CircularAperture([(1028., 1028.)], r=10.)

    # Create the instrument, explicitely set spline order to be used
    cmd = sim.UserCommands(use_instrument=instrument,
                           set_modes=[instrument_mode])
    cmd["!SIM.computing.spline_order"] = spline_order
    opttrain = sim.OpticalTrain(cmd)
    opttrain["detector_linearity"].include = False

    # A progress bar - works only in jupyter notebook
    # pbar = IntProgress(min=0, max=len(cdeltarr))  # instantiate the bar
    # display(pbar)

    # Run the loop
    for cdelt in cdeltarr:
        # Create an image of a point source - a single bright pixel. Header has arbitrary CDELT1/2 to start with.
        img = np.zeros((21, 21), dtype=np.float32)
        img[10, 10] = 1.
        header = fits.Header()
        header.update(
            {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CDELT1': cdelt,
             'CDELT2': cdelt,
             'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec', 'CRPIX1': 10, 'CRPIX2': 10,
             'CRVAL1': 0., 'CRVAL2': 0.,
             'BUNIT': 'Jy'})
        img_hdu = fits.ImageHDU(data=img, header=header)

        src = sim.Source(spectra=[spec], image_hdu=img_hdu)

        opttrain.observe(src, update=True)
        result = opttrain.image_planes[0].data

        starmax.append(result.max())
        bglevel.append(np.median(result[200:800, 200:800]))
        bgnoise.append(np.std(result[200:300, 200:300]))

        phot_table = aperture_photometry(result - bglevel[-1], aperture,
                                         error=np.ones_like(result) * bgnoise[
                                             -1])

        starsum.append(phot_table['aperture_sum'][0])
        starnoise.append(phot_table['aperture_sum_err'][0])

        #pbar.value += 1

    return {"cdelt": cdeltarr, "starsum": starsum, "starnoise": starnoise,
            "bglevel": bglevel,
            "bgnoise": bgnoise, "starmax": starmax}


class TestWhatsGoingOn:
    def test_run_olivers_working_script(self):
        rc.__config__["!SIM.file.local_packages_path"] = r"F:/Work/irdb"
        result = point_source_photometry_versus_cdelt_oliver(3)
        print(result)

    def test_run_miguels_working_script(self):
        rc.__config__["!SIM.file.local_packages_path"] = r"F:/Work/irdb"
        result = point_source_photometry_versus_cdelt_oliver(3)
        print(result)