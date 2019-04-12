import numpy as np
from astropy.io import fits


def _basic_circular_fvpsf():
    psfs = np.zeros([5, 11, 11])
    psfs[0, 1:10, 1:10] = 1
    psfs[1, 2:9, 2:9] = 1
    psfs[2, 3:8, 3:8] = 1
    psfs[3, 4:7, 4:7] = 1
    psfs[4, 5, 5] = 1

    dic = {"CDELT1": 1/3600., "CDELT2":  1/3600., "CRVAL1": 0, "CRVAL2": 0,
           "CRPIX1": 5, "CRPIX2": 5, "CUNIT1": "deg", "CUNIT2": "deg",
           "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "WAVE0": 1.5}
    hdr = fits.Header()
    hdr.update(dic)
    psf_hdu = fits.ImageHDU(data=psfs, header=hdr)

    srmap = np.zeros((100, 100))
    for r in [5, 10, 20, 40]:
        for x in range(100):
            for y in range(100):
                if (x-50)**2 + (y-50)**2 < r**2:
                    srmap[x, y] += 1

    dic = {"CDELT1": 1 / 3600., "CDELT2": 1 / 3600., "CRVAL1": 0, "CRVAL2": 0,
           "CRPIX1": 50, "CRPIX2": 50, "CUNIT1": "deg", "CUNIT2": "deg",
           "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"}
    hdr = fits.Header()
    hdr.update(dic)
    strehl_hdu = fits.ImageHDU(data=srmap, header=hdr)

    dic = {"ECAT": 1, "EDATA": 2}
    hdr = fits.Header()
    hdr.update(dic)
    pri_hdr = fits.PrimaryHDU(header=hdr)

    hdu_list = fits.HDUList([pri_hdr, strehl_hdu, psf_hdu])

    return hdu_list

# psf = _basic_circular_fvpsf()
# psf.writeto("st_circular_fvpsf.fits")
