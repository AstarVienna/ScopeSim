Basic Source from ImageHDUs
===========================

TL;DR
-----

.. jupyter-execute::

    import scopesim, scipy, astropy.io.fits as fits

    hdu = fits.ImageHDU(data=scipy.misc.face(gray=True))
    hdu.header.update({"CDELT1": 1, "CUNIT1": "arcsec", "CRPIX1": 0, "CRVAL1": 0,
                       "CDELT2": 1, "CUNIT2": "arcsec", "CRPIX2": 0, "CRVAL2": 0,})

    ab_spec = scopesim.source.source_templates.ab_spectrum(mag=20)

    image_source = scopesim.Source(image_hdu=hdu, spectra=[ab_spec])

    print(image_source.fields)
    print(image_source.spectra)


Explanation
-----------

.. jupyter-execute::

    import matplotlib.pyplot as plt
    %matplotlib inline

    plt.subplot(121)
    wave = range(3000, 25000)
    plt.plot(wave, image_source.spectra[0](wave))
    plt.xlabel("Wavelength [Angstrom]")
    plt.ylabel("Flux [ph/s/cm2/Angstrom]")
    plt.subplot(122)
    plt.imshow(image_source.fields[0].data)
