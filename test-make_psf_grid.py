import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits

# Create PSF grid
n_psf = 11
s_psf = 128
epsf = [
    [Gaussian2DKernel(1 + 0.1*x + 0.1*y, x_size=s_psf, y_size=s_psf).array
     for y in range(n_psf)] for x in range(n_psf)]
epsf = np.array(epsf)

pix_size = 1.5  # mas
fov = pix_size * 256

grid_spacing = fov/(n_psf-1)
waveleng = 3.7  # um

primary_hdr = fits.Header()
primary_hdr["SIMPLE"] = (True, "conforms to FITS standard")
primary_hdr["BITPIX"] = (8, "array data type")
primary_hdr["NAXIS"] = (0, "number of array dimensions")
primary_hdr["EXTEND"] = True
primary_hdr["FILETYPE"] = 'Point Spread Function (Grid)'
primary_hdr["AUTHOR"] = 'J. Aveiro'
primary_hdr["DATE"] = '2023'
primary_hdr["SOURCE"] = 'TEST'
primary_hdr["ORIGDATE"] = '2023'
primary_hdr["WAVELENG"] = (waveleng, "microns")
primary_hdr["PIXSIZE"] = (pix_size, "milliarcsec")
primary_hdr["XPOSITIO"] = (0.00000, "arcsec")
primary_hdr["YPOSITIO"] = (0.00000, "arcsec")

image_hdr = fits.Header()
image_hdr["WAVELENG"] = (waveleng, "microns")
image_hdr["PIXSIZE"] = (pix_size, "milliarcsec")
image_hdr["NAXIS"] = 4
image_hdr["NAXIS1"] = n_psf
image_hdr["NAXIS2"] = n_psf
image_hdr["NAXIS3"] = s_psf
image_hdr["NAXIS4"] = s_psf
image_hdr["PIXSCALE"] = (pix_size, "milliarcsec")
image_hdr["CDELT1"] = (grid_spacing, "[mas] Coordinate increment at reference point")
image_hdr["CDELT2"] = (grid_spacing, "[mas] Coordinate increment at reference point")
image_hdr["CDELT3"] = (pix_size, "[mas] Coordinate increment at reference point")
image_hdr["CDELT4"] = (pix_size, "[mas] Coordinate increment at reference point")
image_hdr["CTYPE1"] = ("LINEAR", "Coordinate type code")
image_hdr["CTYPE2"] = ("LINEAR", "Coordinate type code")
image_hdr["CTYPE3"] = ("LINEAR", "Coordinate type code")
image_hdr["CTYPE4"] = ("LINEAR", "Coordinate type code")
image_hdr["CUNIT1"] = ("mas", "Units of coordinate increment and value")
image_hdr["CUNIT2"] = ("mas", "Units of coordinate increment and value")
image_hdr["CUNIT3"] = ("mas", "Units of coordinate increment and value")
image_hdr["CUNIT4"] = ("mas", "Units of coordinate increment and value")
image_hdr["CRVAL1"] = (0.0, "[mas] Coordinate value at reference point")
image_hdr["CRVAL2"] = (0.0, "[mas] Coordinate value at reference point")
image_hdr["CRVAL3"] = (0.0, "[mas] Coordinate value at reference point")
image_hdr["CRVAL4"] = (0.0, "[mas] Coordinate value at reference point")
image_hdr["CRPIX1"] = (n_psf//2, "Grid coordinate of reference point")
image_hdr["CRPIX2"] = (n_psf//2, "Grid coordinate of reference point")
image_hdr["CRPIX3"] = (s_psf//2, "Pixel coordinate of reference point")
image_hdr["CRPIX4"] = (s_psf//2, "Pixel coordinate of reference point")

# Construct FITS
primary_hdu = fits.PrimaryHDU(header=primary_hdr)
image_hdu = fits.ImageHDU(epsf, header=image_hdr)
hdul = fits.HDUList([primary_hdu, image_hdu])

# Save
filename = "psf_grid.fits"
hdul.writeto(filename, overwrite=True)
