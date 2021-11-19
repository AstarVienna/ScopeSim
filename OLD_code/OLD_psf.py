"""
PSFs and PSFCubes
=================

.. todo::
    revise this opening text

Description
-----------

.. admonition:: Car Sagan said

    "If you want to bake an apple pie from scratch,
        first you must create the universe"

Single PSFs
-----------

We need to start by generating a single PSF in order to generate a PSFCube.
We need to know the spatial characteristics of the PSF:
The commonalities of all PSFs are:

- pix_width
- pix_height
- pix_res
- type

The types of PSF offered: Moffat, Gaussian2D, Airy, Delta, Line, User
For each of the PSF types we need to create a subclass of PSF. Each subclass
takes its own list of parameters:

- MoffatPSF      (alpha, beta)
- GaussianPSF    (fwhm, eccentricity=0, angle=0)
- AiryPSF        (first_zero, eccentricity=0, angle=0)
- DeltaPSF       (x=0, y=0)
- LinePSF        (x0, x1, y0, y1, angle=0)
- UserPSFCube        (filename, ext_no=0)


Multiple PSFs in a Cube
-----------------------

To generate a PSF cube we need to know the spectral bins and the type of PSF.
The bins are defined by a central wavelength, however a cube should also
contain the edges of each bin so that transmission and emission can be
re-binned properly.
- lam_bin_centers
- lam_bin_edges
- lam_res


A PSF instance will have these additional arguments:
- array ... a 2D array to hold the PSF image

A psf instance will have these additional arguments:
- cube ... a (l,x,y) 3D array to hold the PSF cube

As far as input goes, psf should be able to accept a dictionary with the
keywords necessary to build the cube.

Notes
-----
All wavelength values are given in [um]
All pixel dimensions are given in [arcsec]
All angles are given in [deg]


Classes
-------
PSF(object)
psf(object)


Subclasses
----------
MoffatPSF(PSF)
GaussianPSF(PSF)
AiryPSF(PSF)
DeltaPSF(PSF)
LinePSF(PSF)
UserPSFCube(PSF)

Deltapsf(psf)
Airypsf(psf)
Gaussianpsf(psf)
Moffatpsf(psf)
CombinedPSFCube(psf)
UserPSFCube(psf)
ADC_psf(psf)


There are two types of psf object here:
- a cube
- a single psf image

The cube is essentially a list of psf images, homogenized in size
Should we have separate classes for these?

Both PSF and psf can be created from a single model or through
convolution of a list of PSF components

"""

import logging
from copy import deepcopy

import numpy as np
import scipy.ndimage.interpolation as spi
from scipy.signal import fftconvolve

from astropy.io import fits
#from astropy import units as u   ## unused (OC)
from astropy.convolution import Moffat2DKernel, Gaussian2DKernel
from astropy.convolution import Kernel2D
from astropy.modeling.core import Fittable2DModel
from astropy.modeling.parameters import Parameter

import scopesim.effects.effects_utils
import scopesim.effects.shifts
from scopesim import utils


## TODO
# - Add a ellipticity to the GaussianPSF
# - Make sure MoffatPSF works

#__all__ = []
__all__ = ["PSF", "PSFCube",
           "MoffatPSF", "MoffatPSFCube",
           "AiryPSF", "AiryPSFCube",
           "GaussianPSF", "GaussianPSFCube",
           "DeltaPSF", "DeltaPSFCube",
           "CombinedPSF", "CombinedPSFCube",
           "UserPSF", "UserPSFCube",
           #"poppy_eelt_psf", "poppy_ao_psf", "seeing_psf",
           "get_eelt_segments", "make_foreign_PSF_cube"
           ]



###############################################################################
#                            PSF and PSF subclasses                           #
###############################################################################

# (Sub)Class    Needed              Optional
# PSF           size, pix_res
# DeltaPSF                          pix_res=0.004, size=f(position), position=(0,0)
# AiryPSF       fwhm                pix_res=0.004, size=f(fwhm)
# GaussianPSF   fwhm                pix_res=0.004, size=f(fwhm)
# MoffatPSF     fwhm                pix_res=0.004, size=f(fwhm)
# CombinedPSFCube   psf_list            size=f(psf_list)
# UserPSFCube       filename,           pix_res=0.004, size=f(filename), fits_ext=0


class PSF(object):
    """Point spread function (single layer) base class

    Parameters
    ----------
    size : int
        [pixel] the side length of the array
    pix_res : float
        [arcsec] the pixel scale of the array
    """

    def __init__(self, size, pix_res):
        self.size = size
        self.shape = [size, size]
        self.pix_res = pix_res
        self.array = np.zeros((self.size, self.size))

        self.info = dict([])
        self.info['description'] = "Point spread function (single layer)"

    def __str__(self):
        return self.info['description']

    def set_array(self, array, threshold=1e-15):
        """
        Set the spatial flux distribution array for the PSF

        Renormalise the array make sure there aren't any negative values that
        will screw up the flux

        Parameters
        ----------
        array : np.ndarray
            the array representing the PSF
        threshold : float
            by default set to 1E-15. Below this, the array is set to 0

        """
        self.array = array.astype(np.float32)
        self.array[self.array <= 0] = threshold
        self.array = self.array / np.sum(self.array)
        self.size = self.array.shape[0]
        self.shape = self.array.shape

    def resize(self, new_size):
        """
        Resize the PSF. The target shape is (new_size, new_size).

        Parameters
        ----------
        new_size : int
            the new square dimensions of the PSF array in pixels
        """

         # make sure the new size is always an odd number
        if new_size % 2 == 0:
            new_size += 1

        arr_tmp = np.zeros((new_size, new_size))
        arr_tmp[new_size//2, new_size//2] = 1
        self.set_array(fftconvolve(arr_tmp, self.array, mode="same"))
        #self.set_array(convolve_fft(arr_tmp, self.array))

    def resample(self, new_pix_res):
        """
        Resample the PSF array onto a new grid

        Not perfect, but conserves flux

        Parameters
        ----------
        new_pix_res : float
            [arcsec] the pixel resolution of the returned array

        Returns
        -------
        psf_new : PSF

        Examples
        --------

            >>> new_PSF = old_PSF.resample(new_pix_res)

        """
        scale_factor = self.pix_res / np.float(new_pix_res)
        new_arr = spi.zoom(self.array, scale_factor, order=1)
        new_arr *= np.sum(self.array) / np.sum(new_arr)

        # catch any even-numbered arrays
        if new_arr.shape[0] % 2 == 0:
            tmp_arr = np.zeros((new_arr.shape[0]+1, new_arr.shape[1]+1))
            tmp_arr[:new_arr.shape[0], :new_arr.shape[1]] = new_arr
            new_arr = spi.shift(tmp_arr, 0.5, order=5)

        ############################################################
        # Not happy with the way the returned type is not the same #
        # as the original type. The new object is a plain PSF      #
        ############################################################
        psf_new = PSF(size=new_arr.shape[0], pix_res=new_pix_res)
        psf_new.set_array(new_arr)
        return psf_new

    def convolve(self, kernel):
        """
        Convolve the PSF with another kernel. The PSF keeps its shape

        Parameters
        ----------
        kernel : np.array, PSF
            Either a numpy.ndarray or a PSF (sub)class

        Returns
        -------
        psf_new : PSF

        """
        psf_new = deepcopy(self)

        if issubclass(type(kernel), PSF):
            psf_new.set_array(fftconvolve(self.array, kernel.array, mode="same"))
            #psf_new.set_array(convolve_fft(self.array, kernel.array))
        else:
            psf_new.set_array(fftconvolve(self.array, kernel, mode="same"))
            #psf_new.set_array(convolve_fft(self.array, kernel))
        psf_new.info["Type"] = "Combined"

        return psf_new

    def __array__(self):
        return self.array

    def __mul__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array * x

    def __add__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array + x

    def __sub__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array - x

    def __rmul__(self, x):
        return self.__mul__(x)

    def __radd__(self, x):
        return self.__add__(x)

    def __rsub__(self, x):
        psf_new = deepcopy(self)
        return x - psf_new.array

    def __imul__(self, x):
        return self.__mul__(x)

    def __iadd__(self, x):
        return self.__add__(x)

    def __isub__(self, x):
        return self.__sub__(x)


class DeltaPSF(PSF):
    """
    Generate a PSF with a delta function at position (x,y)

    Parameters
    ----------
    position : tuple
        [pixel] where (x,y) on the array is where the delta function goes
        default is (x,y) = (0,0) and is the centre of the array
    size : int
        [pixel] the side length of the array
    pix_res : float
        [arcsec] the pixel scale used in the array, default is 0.004

    """

    def __init__(self, **kwargs):

        if "position" in kwargs.keys():
            self.position = kwargs["position"]
        else:
            self.position = (0, 0)

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 5
        else:
            size = int(np.max(np.abs(self.position))) * 2 + 5

        if not np.max(self.position) < size:
            raise ValueError("positions are outside array borders:")

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]
        else:
            pix_res = 0.004

        super(DeltaPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Delta"
        self.info["description"] = "Delta PSF, centred at (%.1f, %.1f)" \
                                    % self.position
        self.info["x_shift"] = self.position[0]
        self.info["y_shift"] = self.position[1]

        self.x = self.size // 2 + self.position[0]
        self.y = self.size // 2 + self.position[1]

        x2 = self.x - int(self.x)
        x1 = 1. - x2
        y2 = self.y - int(self.y)
        y1 = 1. - y2

        arr = np.zeros((self.size, self.size))
        arr[int(self.y) : int(self.y) + 2, int(self.x) : int(self.x) + 2] = \
            np.array([[x1 * y1, x2 * y1], [x1 * y2, x2 * y2]])
        self.set_array(arr)



class AiryPSF(PSF):
    """
    Generate a PSF for an Airy function with an equivalent FWHM

    Parameters
    ----------
    fwhm : float
        [arcsec] the equivalent FWHM of the Airy disk core.
    size : int
        [pixel] the side length of the array
    pix_res : float
        [arcsec] the pixel scale used in the array, default is 0.004
    obscuration : float
        [0..1] radius of inner obscuration as fraction of aperture radius
    mode : str
        ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)

    """

    def __init__(self, fwhm, obscuration=0., size=255, pix_res=0.004,
                 **kwargs):

        # Ensure that PSF size is at least 8 times fwhm
        size = int(np.max((round(8 * fwhm / pix_res) * 2 + 1, size)))

        # if size > 511:
            # size = 511
            # print("FWHM [arcsec]:", fwhm, "- pixel res [arcsec]:", pix_res)
            # print("Array size:", size, "x", size, "- PSF FoV:", size * pix_res)
            # logging.warning("PSF dimensions too large - cropped to 512x512")

        # Check for 'mode' keyword argument
        if "mode" in kwargs.keys():
            mode = kwargs["mode"]
        else:
            if size > 100:
                mode = "linear_interp"
            else:
                mode = 'oversample'

        # Ensure that obscuration is between 0 and 1
        if obscuration < 0 or obscuration > 1:
            print("Warning: Obscuration must be between 0 and 1. Using 0.")
            obscuration = 0.

        super(AiryPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Airy"
        self.info['description'] \
            = "Airy PSF, FWHM = %.1f mas, obscuration = %.2f" \
            % (fwhm * 1E3, obscuration)
        self.info["fwhm"] = fwhm * 1E3           # milliarcseconds
        self.info["obscuration"] = obscuration
        self.info["pixel scale"] = pix_res
        self.info["mode"] = mode

        self.fwhm = fwhm
        self.pix_res = pix_res

        # These are first zero and FWHM of J1(x)/x.
        airy_zero = 3.8317059860286755
        airy_fwhm = 2 * 1.616339948310703

        # rescale radial array so that the PSF has the required FWHM
        # (for eps = 0)
        first_zero = fwhm / airy_fwhm * airy_zero / pix_res
        psf_arr = AiryDiskDiff2DKernel(first_zero,
                                       obscuration=obscuration,
                                       x_size=size,
                                       y_size=size,
                                       mode=mode).array
        self.set_array(psf_arr)



class GaussianPSF(PSF):
    """
    Generate a PSF for an Gaussian function

    Parameters
    ----------
    fwhm : float
        [arcsec] the equivalent FWHM of the Airy disk core.
    size : int
        [pixel] the side length of the array
    pix_res : float
        [arcsec] the pixel scale used in the array, default is 0.004
    mode : str
        ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)

    """

    def __init__(self, fwhm, **kwargs):

        self.fwhm = fwhm

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]
        else:
            pix_res = 0.004

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else:
            size = 1

        if "undersized" in kwargs.keys() and kwargs["undersized"]:
            pass
        else:
            size = int(np.max((round(5 * self.fwhm / pix_res) * 2 + 1, size)))

        # Check for 'mode' keyword argument
        if "mode" in kwargs.keys():
            mode = kwargs["mode"]
        else:
            if size > 100:
                mode = "linear_interp"
            else:
                mode = 'oversample'

        super(GaussianPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Gaussian"
        self.info['description'] = "Gaussian PSF, FWHM = %.1f mas" \
                                    % (self.fwhm * 1E3)
        self.info["fwhm"] = self.fwhm * 1E3

        n = (self.fwhm / 2.35) / self.pix_res
        k = Gaussian2DKernel(n, x_size=self.size, y_size=self.size, mode=mode).array

        self.set_array(k)



class MoffatPSF(PSF):
    """
    Generate a PSF for a Moffat function. Alpha is generated from the FWHM and
    Beta = 4.765 (from Trujillo et al. 2001)

    Parameters
    ----------
    fwhm : float
        [arcsec] the equivalent FWHM of the Airy disk core.
    size : int
        [pixel] the side length of the array
    pix_res : float
        [arcsec] the pixel scale used in the array, default is 0.004
    mode : str
        ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)

    """
    def __init__(self, fwhm, **kwargs):

        self.fwhm = fwhm

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]
        else:
            pix_res = 0.004

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else:
            size = 1
        size = int(np.max((round(4 * self.fwhm / pix_res) * 2 + 1, size)))

        # Check for 'mode' keyword argument
        if "mode" in kwargs.keys():
            mode = kwargs["mode"]
        else:
            if size > 100:
                mode = "linear_interp"
            else:
                mode = 'oversample'

        super(MoffatPSF, self).__init__(size, pix_res)

        beta = 4.765 ### Trujillo et al. 2001
        alpha = self.fwhm/(2 * np.sqrt(2**(1/beta) - 1))
        self.info["Type"] = "Moffat"
        self.info['description'] = "Moffat PSF, FWHM = %.1f, alpha = %.1f"\
                                       % (self.fwhm * 1E3, alpha)
        self.info["fwhm"] = self.fwhm * 1E3

        self.set_array(Moffat2DKernel(alpha, beta, x_size=self.size,
                                      y_size=self.size, mode=mode).array)



class CombinedPSF(PSF):
    """
    Generate a PSF from a collection of several PSFs.

    Parameters
    ----------
    psf_list : list
        A list of PSF py_objects
    size : int
        [pixel] the side length in pixels of the array

    """

    def __init__(self, psf_list, **kwargs):
        """Generate a master psf through convolution of a list of psfs"""

        if not hasattr(psf_list, "__len__") or len(psf_list) < 2:
            raise ValueError("psf_list requires more than 1 PSF object")

        pix_res_list = [psf.pix_res for psf in psf_list]
        if not all(res == pix_res_list[0] for res in pix_res_list):
            raise ValueError("Not all PSFs have the same pixel resolution")

        pix_res = pix_res_list[0]

        if "size" in kwargs.keys():
            size = int(kwargs["size"] // 2) * 2 + 1
        else:
            size_list = [psf.size for psf in psf_list]
            size = int(np.max(size_list) // 2) * 2 + 1

        ## Compensate for the shift in centre due to a DeltaPSF
        shifts = np.asarray([(0, 0)] + [psf.position for psf in psf_list \
                                        if psf.info["Type"] == "Delta"])
        size += 2 * np.max(shifts)

        arr_tmp = np.zeros((size, size))
        arr_tmp[size // 2, size // 2] = 1

        for psf in psf_list:
            arr_tmp = fftconvolve(arr_tmp, psf.array, mode="same")
            #arr_tmp = convolve_fft(arr_tmp, psf.array)


        super(CombinedPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Combined"
        self.info['description'] = "Combined PSF from " + str(len(psf_list)) \
                                                                + "PSF py_objects"
        self.set_array(arr_tmp)



class UserPSF(PSF):
    """
    Import a PSF from a FITS file.

    Parameters
    ----------
    filename : str
        path to the FITS file to be read in
    fits_ext : int, optional
        the FITS extension number (default 0) for the data
    pix_res : float, optional
        [arcsec] the pixel scale used in the array, default is 0.004

    """

    def __init__(self, filename, **kwargs):

        if "fits_ext" in kwargs.keys():
            fits_ext = kwargs["fits_ext"]
        else:
            fits_ext = 0

        self.filename = filename
        self.fits_ext = fits_ext

        header = fits.getheader(self.filename, ext=self.fits_ext)
        data = fits.getdata(self.filename, ext=self.fits_ext)
        size = header["NAXIS1"]

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]
        elif "CDELT1" in header.keys():
            pix_res = header["CDELT1"]
        elif "CD1_1" in header.keys():
            pix_res = header["CD1_1"]
        else:
            pix_res = 0.004

        # If pix_res is
        #   * >1E-1 it must be in mas
        #   * ~1E-3, it must be in arcsec
        #   * <1E-5 it must be in deg
        if pix_res > 1E-1:
            pix_res *= 1E-3
        elif pix_res < 1E-5:
            pix_res *= 3600

        super(UserPSF, self).__init__(size, pix_res)
        self.info["Type"] = "User"
        self.info['description'] = "PSF from FITS file: " + self.filename

        self.set_array(data)

        if "size" in kwargs.keys():
            self.resize(kwargs["size"])



###############################################################################
#                       psf and psf subclasses                        #
###############################################################################

# (Sub)Class    Needed              Optional
# PSF           size, pix_res
# DeltaPSF                          pix_res=0.004, size=f(position), position=(0,0)
# AiryPSF       fwhm                pix_res=0.004, size=f(fwhm)
# GaussianPSF   fwhm                pix_res=0.004, size=f(fwhm)
# MoffatPSF     fwhm                pix_res=0.004, size=f(fwhm)
# CombinedPSFCube   psf_list            size=f(psf_list)
# UserPSFCube       filename,           pix_res=0.004, size=f(filename), fits_ext=0

#    Keywords:
#    - type:
#    - lam_bin_centers

#    Bound keywords:
#    - fwhm : needed for AiryPSF, GaussianPSF, MoffatPSF
#    - psf_list: needed for CombinedPSFCube
#    - filename: needed for UserPSFCube

#    Optional keywords
#    - size:
#    - pix_res:
#    - position: optional in DeltaPSF
#    - fits_ext: optional in UserPSFCube


class PSFCube(object):
    """
    Class holding wavelength dependent point spread function.

    Parameters
    ----------
    lam_bin_centers : array
        [um] the centre of each wavelength slice

    Notes
    -----
    - len(self) return the number of layers in the psf
    - self[i] returns the PSF object for layer i. If the __array__ function
      is called, self[i] will return the array associated with the instance
      e.g plt.imshow(self[i]) will plot PSF.array from self.psf_slices[i]
    - Maths operators \*,\+,\- act equally on all PSF.arrays in self.psf_slices


    """

    def __init__(self, lam_bin_centers):

        self.lam_bin_centers = lam_bin_centers
        self.psf_slices = [None] * len(lam_bin_centers)

        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "Point spread function (multiple layer)"


    def resize(self, new_size):
        """
        Resize the list of PSFs. The target shape is (new_size, new_size).

        Parameters
        ----------
        new_size : int
            [pixel] the new size of the PSF array in pixels

        """
        for psf in self.psf_slices:
            psf.resize(new_size)

    def resample(self, new_pix_res):
        """
        Resample the the list of PSF array onto a new grid

        Not perfect, but conserves flux

        Parameters
        ----------
        new_pix_res : float
            [arcsec] the pixel resolution of the returned array

        Returns
        -------
        psf_new : PSF

        Examples
        --------

            >>> new_PSF = old_PSF.resample(new_pix_res)

        """

        ## TODO: Check whether this makes sense
        self.psf_slices = [psf.resample(new_pix_res) for psf in self.psf_slices]

    def export_to_fits(self, filename, clobber=True):
        """
        Export the psf to a FITS file for later use

        Parameters
        ----------
        filename : str

        """

        ## TODO: use header_info or drop it (OC)
        ## KL: done

        ext_list = fits.HDUList()

        #for i in range(len(self)):
        #    psf = self.psf_slices[i]
        for i, psf in enumerate(self.psf_slices):
            #if i == 0:
            #    hdu = fits.PrimaryHDU(psf.array)
            #else:
            ## PrimaryHDUs are typically empty. Not necessary (OC)
            hdu = fits.ImageHDU(psf.array)
            hdu.header["CDELT1"] = (psf.pix_res, "[arcsec] - Pixel resolution")
            hdu.header["CDELT2"] = (psf.pix_res, "[arcsec] - Pixel resolution")
            hdu.header["WAVE0"] = (self.lam_bin_centers[i],
                                   "[micron] - Wavelength of slice")
            hdu.header["NSLICES"] = (len(self), "Number of wavelength slices")

            for k in psf.info.keys():
                hdu.header[k[:8].upper()] = (self[i].info[k], k)

            ext_list.append(hdu)

        ext_list.writeto(filename, clobber=clobber, checksum=True)


    def convolve(self, kernel_list):
        """
        Convolve a list of PSFs with a list of kernels

        Parameters
        ----------
        kernel_list : list
            list of PSF py_objects of 2D arrays

        """
        if len(self.psf_slices) != len(kernel_list):
            print("len(PSF_slices):", len(self.psf_slices),
                  "len(kernel_list):", len(kernel_list))
            raise ValueError("Number of kernels must equal number of PSFs")

        for i in np.arange(len(self.psf_slices)):
            tmp = self.psf_slices[i].convolve(kernel_list[i])
            self.psf_slices[i].set_array(tmp.array)
        self.info["Type"] = "Complex"


    def nearest(self, lam):
        """
        Returns the PSF closest to the desired wavelength, lam [um]

        Parameters
        ----------
        lam : float
            [um] desired wavelength

        Returns
        -------
        psf_slice : PSF

        """

        i = utils.nearest(self.lam_bin_centers, lam)
        return self.psf_slices[i]


    def __str__(self):
        return self.info['description']

    def __getitem__(self, i):
        if len(self.psf_slices) > 1:
            return self.psf_slices[i]
        else:
            return self.psf_slices[0]

    def __len__(self):
        return len(self.psf_slices)


    def __mul__(self, x):
        newpsf = deepcopy(self)

        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else:
            y = x

        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")

        for i in np.arange(len(self.psf_slices)):
            newpsf[i].set_array(self.psf_slices[i] * y[i])
        return newpsf



    def __add__(self, x):
        newpsf = deepcopy(self)

        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else:
            y = x

        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")

        for i in np.arange(len(self.psf_slices)):
            newpsf[i].set_array(self.psf_slices[i] + y[i])
        return newpsf


    def __sub__(self, x):
        newpsf = deepcopy(self)

        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else:
            y = x

        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")

        for i in np.arange(len(self.psf_slices)):
            newpsf[i].set_array(self.psf_slices[i] - y[i])
        return newpsf



    def __rmul__(self, x):
        self.__mul__(x)

    def __radd__(self, x):
        self.__add__(x)

    def __rsub__(self, x):
        self.__sub__(x)



class DeltaPSFCube(PSFCube):
    """
    Generate a list of DeltaPSFs for wavelengths defined in lam_bin_centers

    Parameters
    ----------
    lam_bin_centers : array
        [um] the centre of each wavelength slice
    positions : tuple, list
        (x,y) either a tuple, or a list of tuples denoting the
                 position of the delta function
    pix_res : float
        [arcsec], pixel scale of the PSF. Default is 0.004 arcsec
    size : int
        [pixel] side length of the PSF array

    """

    def __init__(self, lam_bin_centers, positions=(0, 0), **kwargs):
        super(DeltaPSFCube, self).__init__(lam_bin_centers)

        if not hasattr(positions[0], "__len__"):
            positions = [positions]*len(self)

        for i in range(len(self)):
            self.psf_slices[i] = DeltaPSF(position=positions[i], **kwargs)

        self.info['description'] = "List of Delta function PSFs"
        self.info["Type"] = "DeltaCube"


class AiryPSFCube(PSFCube):
    """
    Generate a list of AiryPSFs for wavelengths defined in lam_bin_centers

    Parameters
    ----------
    lam_bin_centers : array, list
        [um] a list with the centres of each wavelength slice
    fwhm : array, list
        [arcsec] the equivalent FWHM of the PSF.
    diameter : float
        [m] diamter of primary mirror. Default is 39.3m.
    mode : str
        ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)

    """

    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(AiryPSFCube, self).__init__(lam_bin_centers)

        if "diameter" in kwargs.keys():
            self.diameter = kwargs["diameter"]
        else:
            self.diameter = 39.3

        if "obscuration" in kwargs.keys():
            self.obscuration = kwargs["obscuration"]
        else:
            self.obscuration = 11.1/39.3

        if fwhm is None:
            # lam in um, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)
        else:
            self.fwhm = fwhm

        self.psf_slices = [AiryPSF(fwhm=f, **kwargs) for f in self.fwhm]
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "List of Airy function PSFs"
        self.info["Type"] = "AiryCube"


class GaussianPSFCube(PSFCube):
    """
    Generate a list of GaussianPSFs for wavelengths defined in lam_bin_centers

    Parameters
    ----------
    lam_bin_centers : array, list
        [um] a list with the centres of each wavelength slice
    fwhm : array, list
        [arcsec] the equivalent FWHM of the PSF.
    diameter : float
        [m] diamter of primary mirror. Default is 39.3m.
    mode : str
        ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)

    """

    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(GaussianPSFCube, self).__init__(lam_bin_centers)

        if "diameter" in kwargs.keys():
            self.diameter = kwargs["diameter"]
        else:
            self.diameter = 39.3

        if fwhm is None:
            # lam in um, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)

        self.psf_slices = [GaussianPSF(fwhm=f, **kwargs) for f in self.fwhm]
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "List of Gaussian function PSFs"
        self.info["Type"] = "GaussianCube"

class MoffatPSFCube(PSFCube):
    """
    Generate a list of MoffatPSFs for wavelengths defined in lam_bin_centers

    Parameters
    ----------
    lam_bin_centers : array, list
        [um] a list with the centres of each wavelength slice
    fwhm : array, list
        [arcsec] the equivalent FWHM of the PSF.
    diameter : float
        [m] diamter of primary mirror. Default is 39.3m.
    mode : str
        ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)

    """

    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(MoffatPSFCube, self).__init__(lam_bin_centers)

        if "diameter" in kwargs.keys():
            self.diameter = kwargs["diameter"]
        else:
            self.diameter = 39.3

        if fwhm is None:
            # lam in um, diameter in m, 206265 is 1 rad in arcsec
            rad2arcsec = 3600 * 180. / np.pi
            self.fwhm = rad2arcsec * 1.22 * lam_bin_centers * 1E-6 / self.diameter
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)

        self.psf_slices = [MoffatPSF(fwhm=f, **kwargs) for f in fwhm]
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "List of Moffat function PSFs"
        self.info["Type"] = "MoffatCube"

class CombinedPSFCube(PSFCube):
    """
    Generate a list of CombinedPSFCubes from the list of psfs in psf_list

    Parameters
    ----------
    lam_bin_centers : array, list
        [um] a list with the centres of each wavelength slice
    fwhm : array, list
        [arcsec] the equivalent FWHM of the PSF.
    diameter : float
        [m] diamter of primary mirror. Default is 39.3m.

    """

    def __init__(self, psf_list, **kwargs):

        if not (isinstance(psf_list, list) and len(psf_list) >= 2):
            raise ValueError("psf_list only takes a list of psf py_objects")

        ## Check that the wavelengths are equal
        lam_list = [cube.lam_bin_centers for cube in psf_list]
        if not all([all(lam == lam_list[0]) for lam in lam_list]):
            raise ValueError("Wavelength arrays of psf cubes are not equal")
        lam_bin_centers = lam_list[0]

        super(CombinedPSFCube, self).__init__(lam_bin_centers)

        self.info['description'] = "Master psf cube from list"
        self.info["Type"] = "CombinedCube"

        for i, psfi in enumerate(psf_list):
            self.info['PSF%02d' % (i+1)] = psfi.info['description']

        for i in range(len(self)):
            self.psf_slices[i] = CombinedPSFCube([psf[i] for psf in psf_list],
                                                 **kwargs)
        self.size = [psf.size for psf in self.psf_slices]


class UserPSFCube(PSFCube):
    """
    Read in a psf previously saved as a FITS file

    Keywords needed for a psf to be read in:
    NSLICES, WAVECENT, NAXIS1, CDELT1, PSF_TYPE, DESCRIPT

    Parameters
    ----------
    filename : str
        the path to the FITS file holding the cube

    Notes
    -----
    A separate function will exist to convert foreign PSF FITS files into
    ``scopesim.psf`` readable FITS files

    See Also
    --------
    :func:`.make_foreign_PSF_cube`

    """

    def __init__(self, filename, lam_bin_centers):
        if not hasattr(lam_bin_centers, "__len__"):
            lam_bin_centers = [lam_bin_centers]

        psf_slices = []

        # pull out the wavelengths in the PSF FITS files
        psf_lam_cen = []

        if isinstance(filename, fits.HDUList):
            hdulist = filename
        else:
            hdulist = fits.open(filename)

        n_slices = len(hdulist.info(output=False))


        for i in range(n_slices):
            if "WAVE0" in hdulist[i].header.keys():
                psf_lam_cen += [hdulist[i].header["WAVE0"]]
            elif "WAVELENG" in hdulist[i].header.keys():
                psf_lam_cen += [hdulist[i].header["WAVELENG"]]
            else:
                raise ValueError("""Could not determine wavelength of PSF in
                                 extension """ + str(i) + """. FITS file
                                 needs either WAVE0 or WAVELENG header
                                 keywords. \n Use scopesim.utils.add_keyword()
                                 to add WAVELENG to the FITS header""")
            # If the wavelength is not in the 0.1-2.5 range, it must be in [m]
            if psf_lam_cen[i] < 0.1:
                psf_lam_cen[i] *= 1E6


        # find the closest PSFs in the file to what is needed for the psf
        i_slices = utils.nearest(psf_lam_cen, lam_bin_centers)
        i_psf_lam_cen = [psf_lam_cen[i] for i in np.unique(i_slices)]

        # import only the relevant PSFs
        for i in np.unique(i_slices):

            hdr = hdulist[i].header
            self.header = hdr

            if 'CDELT1' in hdr.keys():
                pix_res = hdr["CDELT1"]
            elif 'CD1_1' in hdr.keys():
                pix_res = hdr['CD1_1']
            elif 'PIXSCALE' in hdr.keys():
                pix_res = hdr['PIXSCALE']
            else:
                raise KeyError("Could not get pixel scale from " +
                               filename)

            if pix_res > 1:
                logging.warning("CDELT > 1. Assuming the scale to be [mas]")
                pix_res *= 1E-3

            psf = PSF(size=hdr["NAXIS1"], pix_res=pix_res)
            psf.set_array(hdulist[i].data)

            if "PSF_TYPE" in hdr.keys():
                psf.info["Type"] = hdr["PSF_TYPE"]
            else:
                psf.info["Type"] = "Unknown"

            if "DESCRIPT" in hdr.keys():
                psf.info["description"] = hdr["DESCRIPT"]
            else:
                psf.info["description"] = "Unknown"

            psf_slices += [psf]


        hdulist.close()

        lam_bin_centers = np.array(lam_bin_centers)

        super(UserPSFCube, self).__init__(i_psf_lam_cen)
        self.psf_slices = psf_slices
        self.size = [psf.size for psf in self.psf_slices]

        if isinstance(filename, str):
            self.info['description'] = "User PSF cube input from " + filename
        else:
            self.info['description'] = "User PSF cube input from memory"
        self.info["Type"] = psf_slices[0].info["Type"]+"Cube"



class ADC_PSFCube(DeltaPSFCube):
    """
    Generates a DeltaPSFCube with the shifts required to mimic the ADC

    Parameters
    lam_bin_centers : array, list
        [um] a list with the centres of each wavelength slice

    Other Parameters:
    pix_res : flaot
        [arcsec] the pixel scale used in the array, default is 0.004
    OBS_PARALLACTIC_ANGLE : float
        [deg] the orientation of the input cube relative to the zenith
    INST_ADC_EFFICIENCY : float
        [%] efficiency of the ADC
    SCOPE_LATITUDE : float
        [deg] latitude of the telescope site in decimal degrees
    SCOPE_ALTITUDE : float
        [m] height above sea level of the telescope site
    ATMO_REL_HUMIDITY : float
        [%] relative humidity in percent
    ATMO_AIRMASS : float
        airmass of the object
    ATMO_TEMPERATURE : float
        [deg C] air temperature of the observing site in Celsius
    ATMO_PRESSURE : float
        [mbar] air pressure of the observing site in millibar

    The default values for the above mentioned keywords are:
    0.004 arcsec, 0 deg, 100%, -24.5 deg, 3064m, 60%, 60 deg, 0 C, 750mbar
    """

    def __init__(self, lam_bin_centers, **kwargs):
        params = {"pix_res"           :0.004,
                  "PARALLACTIC_ANGLE"  :0,
                  "INST_ADC_EFFICIENCY":100,
                  "SCOPE_LATITUDE"     :-24.5,
                  "SCOPE_ALTITUDE"     :3064,
                  "ATMO_REL_HUMIDITY"  :60,
                  "ATMO_AIRMASS"       :2.,
                  "ATMO_TEMPERATURE"   :0,
                  "ATMO_PRESSURE"      :750}

        params.update(**kwargs)
        pix_res = params["pix_res"]
        para_angle = params["OBS_PARALLACTIC_ANGLE"]
        effectiveness = params["INST_ADC_PERFORMANCE"] / 100.

        ## get the angle shift for each slice
        zenith_distance = utils.airmass2zendist(params["ATMO_AIRMASS"])
        angle_shift = [
            scopesim.effects.shifts.atmospheric_refraction(lam,
                                                           zenith_distance,
                                                           params["ATMO_TEMPERATURE"],
                                                           params["ATMO_REL_HUMIDITY"],
                                                           params["ATMO_PRESSURE"],
                                                           params["SCOPE_LATITUDE"],
                                                           params["SCOPE_ALTITUDE"])
            for lam in lam_bin_centers]

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        pixel_shift = (angle_shift - angle_shift[-1]) / pix_res
        if np.max(np.abs(pixel_shift)) > 1000:
            raise ValueError("Pixel shifts too great (>1000), check units")

        ## Rotate by the paralytic angle
        x = -pixel_shift * np.sin(np.deg2rad(para_angle)) * (1. - effectiveness)
        y = -pixel_shift * np.cos(np.deg2rad(para_angle)) * (1. - effectiveness)
        positions = [(xi, yi) for xi, yi in zip(x, y)]

        super(ADC_PSFCube, self).__init__(lam_bin_centers,
                                          positions=positions,
                                          pix_res=pix_res)
        self.info["Type"] = "ADC_psf"
        self.info['description'] = "ADC PSF cube for ADC effectiveness:" + \
                                    str(params["INST_ADC_EFFICIENCY"]) + \
                                    ", z0:" + str(params["ATMO_AIRMASS"])


# The following two classes implement a kernel for the PSF of a centrally
# obscured circular aperture. The classes are modelled after the kernels
# in astropy.convolution.kernel and the models in astropy.modeling.models,
# both from astropy version 1.1.1.
class AiryDiskDiff2DKernel(Kernel2D):
    """
    2D kernel for PSF for annular aperture

    This kernel models the diffraction pattern of a circular aperture with
    a central circular obscuration. This kernel is normalized to a peak
    value of 1.

    Parameters
    ----------
    radius : float
        The radius of the unobscured Airy disk kernel (radius of the first
        zero). Compute this from the outer aperture radius.
    obscuration : float
        Fraction of the aperture that is obscured (inner radius / outer radius)
        Default obscuration = 0.
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = 8 * radius.
    y_size : odd int, optional
        Size in y direction of the kernel array. Default = 8 * radius.
    mode : str, optional
        One of the following discretization modes:
            * 'center' (default)
                Discretize model by taking the value
                at the center of the bin.
            * 'linear_interp'
                Discretize model by performing a bilinear interpolation
                between the values at the corners of the bin.
            * 'oversample'
                Discretize model by taking the average
                on an oversampled grid.
            * 'integrate'
                Discretize model by integrating the
                model over the bin.
    factor : number, optional
        Factor of oversampling. Default factor = 10.

    See Also
    --------
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, MexicanHat2DKernel,
    Ring2DKernel, TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel
    (in astropy.kernels)

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from scopesim.psf import AiryDiskDiff2DKernel
        airydiskdiff_2D_kernel = AiryDiskDiff2DKernel(10)
        plt.imshow(airydiskdiff_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()

    """
    _is_bool = False


    def __init__(self, radius, obscuration, **kwargs):
        from astropy.convolution.kernels import _round_up_to_odd_integer

        self._model = AiryDiskDiff2D(1, 0, 0, radius, obscuration)
        self._default_size = _round_up_to_odd_integer(8 * radius)
        super(AiryDiskDiff2DKernel, self).__init__(**kwargs)
        self.normalize()
        self._truncation = None

class AiryDiskDiff2D(Fittable2DModel):
    """
    Two dimensional Airy disk model with central obscuration.

    Parameters
    ----------
    amplitude : float
        Amplitude of the function.
    x_0 : float
        x position of the maximum of the function
    y_0 : float
        y position of the maximum of the function
    radius : float
        The radius of the unobscured Airy disk (radius of the first zero).
    eps : float
        The ratio of the inner to the outer radius of an annular aperture.

    See Also
    --------
    AiryDisk2D

    Notes
    -----
    Model formula:

        .. math:: f(r) = A \\left[\\frac{2 J_1(\\frac{\\pi r}{R/R_z})}\\right]^2

    Where :math: `J_1` is the first order Bessel function of the first
    kind, :math: `r` is the radial distance from the maximum of the
    function (:math: `r = \\sqrt{(x - x_0)^2 + (y - y_0)^2}`), :math:`R`
    is the input ``radius`` parameter, and :math:`R_z =
    1.2196698912665045`.

    For an optical system, the radius of the first zero represents the
    limiting angular resolution and is approximately 1.22 * lambda / D,
    where lambda is the wavelength of the light and D is the diameter of
    the aperture.
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    radius = Parameter(default=1)
    eps = Parameter(default=0)
    _j1 = None

    def __init__(self, amplitude=amplitude.default, x_0=x_0.default,
                 y_0=y_0.default, radius=radius.default, eps=eps.default,
                 **kwargs):
        if self._j1 is None:
            try:
                from scipy.special import j1, jn_zeros
                self.__class__._j1 = j1
                self.__class__._rz = jn_zeros(1, 1)[0] / np.pi
            #add a ValueError here for python3 + scipy < 0.12
            except ImportError:
                raise ImportError("AiryDiskDiff2D model requires scipy > 0.11.")

        super(AiryDiskDiff2D, self).__init__(amplitude=amplitude,
                                             x_0=x_0,
                                             y_0=y_0,
                                             radius=radius,
                                             eps=eps,
                                             **kwargs)

    # Comment and methods copied from astropy v1.1.1
    # TODO: Why does this particular model have its own special __deepcopy__
    # and __copy__?  If it has anything to do with the use of the j_1 function
    # that should be reworked.
    def __deepcopy__(self, memo):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    def __copy__(self):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    @classmethod
    def evaluate(cls, x, y, amplitude, x_0, y_0, radius, eps):
        """Two dimensional Airy difference model function"""

        r = np.sqrt((x - x_0)**2 + (y - y_0)**2) / (radius / cls._rz)
        # Since r can be zero, we have to take care to treat that case
        # separately so as not to raise a numpy warning
        z = np.ones(r.shape)
        rt = np.pi * r[r > 0]
        z[r > 0] = (2.0 * cls._j1(rt) / rt -
                    2.0 * eps * cls._j1(eps * rt) / rt)**2
        z *= amplitude
        return z


################################################################################
# Convenience functions


def make_foreign_PSF_cube(fnames, out_name=None, window=None, pix_res_orig=None,
                          pix_res_final=None, wavelengths=None):
    """
    Combine several PSF FITS images into a single PSF FITS file

    Parameters
    ----------
    fnames : list
        List of path names to the FITS files
    out_name : str, optional
        If out_name is not ``None``, the resulting FITS file is saved under the
        name ``out_name``
    window : int, list, tuple, optional
        If window is not ``None``, a windowed section of the PSFs are extracted
        window = (left, right, top, bottom)
        window = square radius


    Examples
    --------
    ::

        >>> from glob import glob
        >>> import scopesim as sim
        >>> fnames = glob("D:\Share_VW\Data_for_SimCADO\PSFs\yann_2016_11_10\*.fits")
        >>> sim.psf.make_foreign_PSF_cube(fnames, "PSF_SCAO.fits",
                                          window=512,
                                          pix_res_orig=[0.0028, 0.0037, 0.00492],
                                          pix_res_final=[0.004, 0.004, 0.004],
                                          wavelengths=[1.25,1.65,2.2])

    """

    if pix_res_orig is None:
        pix_res_orig = [None]*len(fnames)

    if pix_res_final is None:
        pix_res_final = [None]*len(fnames)

    if wavelengths is None:
        wavelengths = [None]*len(fnames)


    hdu_list = fits.HDUList()
    for fname, res_i, res_f, wave in zip(fnames, pix_res_orig,
                                         pix_res_final, wavelengths):
        dat = fits.getdata(fname)
        hdr = fits.getheader(fname)

        if res_i is not None and res_f is not None:
            if res_i != res_f:
                dat = spi.zoom(dat, res_i/res_f)

        if window is not None:
            if type(window) in (tuple, list):
                w = window
                dat = dat[w[0]:w[1],w[2]:w[3]]
            elif type(window) == int:
                dw = window
                xc, yc = dat.shape[0] // 2, dat.shape[1] // 2
                dat = dat[xc-dw:xc+dw, yc-dw:yc+dw]
            else:
                print("Unknown type for window", type(window))

        scale_factor = 1./np.sum(dat)
        dat *= scale_factor

        hdu = fits.ImageHDU(dat, hdr)

        if res_f is not None:
            hdu.header["CDELT1"] = (res_f, "[arcsec] pixel resolution")
            hdu.header["CDELT2"] = (res_f, "[arcsec] pixel resolution")
        if wave is not None:
            hdu.header["WAVE0"] = (wave, "[micron] - Wavelength of slice")

        hdu_list.append(hdu)

    if out_name is None:
        return hdu_list
    else:
        hdu_list.writeto(out_name, clobber=True)



def poppy_ao_psf(strehl, mode="wide", plan="A", size=1024, filename=None,
                 **kwargs):
    """
    Create a diffraction limited E-ELT PSF with a Seeing halo

    Uses POPPY to create a diffraction limited PSF for the E-ELT for a
    certain configuration of mirror segments. The diffraction limited core is
    added to seeing halo, modelled by either a moffat or gassian profile.


    Parameters
    ----------
    strehl : float
        [0 .. 1] The components are summed and weighted according to the strehl
        ratio
        psf = (1-strehl)*seeing_psf + (strehl)*diff_limited_psf
    mode : str, optional
        ["wide", "zoom"] Default = "wide". Sets the pixel size for each
        of the MICADO imaging modes - {"wide" : 4mas, "zoom" : 1.5mas}
    plan : str, optional
        ["A", "B"], Default = "A"
        * Plan A is for a fully populated mirror (798 segments)
        * Plan B has the inner 5 rings missing (588 segments) and a further
        5 random segments missing (583 segments)
    size : int, optional
        [pixels] Default = 1024
    filename : str, optional
        Default = None. If filename is not None, the resulting FITS object will
        be saved to disk


    Other Parameters
    ----------------
    fwhm : float
        [arcsec] Default : 0.8
    psf_type : str
        Default : "moffat"
    wavelength : float
        [um] Default  : 2.2
    segments : list
        Default : None. A list of which hexagonal poppy segments to use
        See :func:`.get_eelt_segments`
    flattoflat : float
        [m] Default : 1.256
    gap  : float
        [m] Default  : 0.004
    secondary_radius  : float
        [m] Default : 5
    n_supports : int
        Default : 6
    support_width : float
        [m] Default : 0.2
    support_angle_offset  : float
        [deg] Default : 0
    n_missing : int
        Default : None. Number of segments missing
    pupil_inner_radius  : float
        [m] Default : None  # Plan A: 5.6m, Plan B: 11.5m
    pupil_outer_radius  : float
        [m] Default : 19


    Returns
    -------
    hdu_list : astropy.HDUList
        an astropy FITS object containing the PSFs for the given wavelengths

    See Also
    --------
    :func:`.get_eelt_segments`

    """
    try:
        import poppy
    except:
        raise ImportError("""Poppy is not installed -
        See https://pythonhosted.org/poppy""")

    params = {"strehl"               : strehl,
              "mode"                 : mode,
              "size"                 : size,
              "plan"                 : plan,
              "fwhm"                 : 0.8,
              "psf_type"             : "moffat",
              "wavelength"           : 2.2,
              "segments"             : None,
              "flattoflat"           : 1.256,
              "gap"                  : 0.004,
              "secondary_radius"     : 5,
              "n_supports"           : 6,
              "support_width"        : 0.2,
              "support_angle_offset" : 0,
              "n_missing"            : None,
              "use_pupil_mask"        : True,
              "pupil_inner_radius"    : None,
              "pupil_outer_radius"    : 19
              }
    params.update(kwargs)

    if "pix_res" not in params:
        params["pix_res"] = 0.0015 if mode.lower() == "zoom" else 0.004

    #work out the multiple wavelength
    wavelength = params["wavelength"]
    if np.isscalar(wavelength):
        wavelength = np.array([wavelength])

    # Make the PSF(s) for the E-ELT
    eelt = fits.HDUList()
    for wave in wavelength:
        print("Generating an E-ELT PSF at", wave, "[um]")
        params["wavelength"] = wave
        eelt.append(poppy_eelt_psf(**params)[0])

    # Make a seeing limited PSF
    seeing = seeing_psf(fwhm     = params["fwhm"],
                        psf_type = params["psf_type"],
                        size     = size,
                        pix_res  = params["pix_res"])

    # Combine the two PSF FITS py_objects
    hdu_list = fits.HDUList()
    for psf in eelt:
        poppy_ao = fits.ImageHDU()
        poppy_ao.data = (1. - strehl) * seeing[0].data + strehl * psf.data
        poppy_ao.data = poppy_ao.data.astype(np.float32)

        poppy_ao.header["PIXELSCL"] = params["pix_res"]
        poppy_ao.header["PSF_TYPE"] = "POPPY"
        poppy_ao.header["CDELT1"]   = params["pix_res"]
        poppy_ao.header["CDELT2"]   = params["pix_res"]
        poppy_ao.header["CUNIT1"]   = "arcsec"
        poppy_ao.header["CUNIT2"]   = "arcsec"

        for key in params:
            try: poppy_ao.header[key] = params[key]
            except: pass

        poppy_ao.header.extend(psf.header.cards)

        hdu_list.append(poppy_ao)

    if filename is None:
        return hdu_list
    else:
        print("Writing to", filename)
        hdu_list.writeto(filename, clobber=True)


def seeing_psf(fwhm=0.8, psf_type="moffat", size=1024, pix_res=0.004,
               filename=None):
    """
    Return a seeing limited PSF

    Parameters
    ----------
    fwhm : float, optional
        [arcsec] Default = 0.8
    psf_type : str, optional
        ["moffat, "gaussian"] Default = "moffat"
    size : int, optional
        [pixel] Default = 1024
    pix_res : float, optional
        [arcsec] Default = 0.004
    filename : str, optional
        Default = None. If filename is not None, the resulting FITS object will
        be saved to disk

    Returns
    -------
    seeing_psf : 2D-array

    Notes
    -----
    # Moffat description
    # https://www.gnu.org/software/gnuastro/manual/html_node/PSF.html
    #
    # Approximate parameters - Bendinelli 1988
    # beta = 4.765 - Trujillo et al. 2001

    """

    if fwhm > 5:
        logging.warning("FWHM is rather large: [arcsec]"+str(fwhm))
    fwhm_pix = fwhm/pix_res

    if "moff" in psf_type.lower():
        beta = 4.785
        alpha = fwhm_pix / (2 * np.sqrt(2**(1/beta)-1))

        # astropy gamma = Trujillo alpha
        # astropy alpha = Trujillo beta
        seeing_psf = Moffat2DKernel(gamma=alpha, alpha=beta,
                                    x_size=size, y_size=size,
                                    factor=1).array
    elif "gauss" in psf_type.lower():
        sigma = fwhm_pix/2.3548
        seeing_psf = Gaussian2DKernel(stddev=sigma,
                                      x_size=size, y_size=size,
                                      factor=1).array

    hdu = fits.PrimaryHDU(seeing_psf)
    hdu.header["PIXELSCL"] = pix_res
    hdu.header["PSF_TYPE"] = psf_type
    hdu.header["FWHM"]     = fwhm
    hdu.header["CDELT1"]    = pix_res
    hdu.header["CDELT2"]    = pix_res
    hdu.header["CUNIT1"]    = "arcsec"
    hdu.header["CUNIT2"]    = "arcsec"

    hdu_list = fits.HDUList([hdu])

    if filename is None:
        return hdu_list
    else:
        print("Writing to", filename)
        hdu_list.writeto(filename, clobber=True)


def poppy_eelt_psf(plan="A", wavelength=2.2, mode="wide", size=1024,
                  segments=None, filename=None, use_pupil_mask=True, **kwargs):
    """
    Generate a PSF for the E-ELT for plan A or B with POPPY

    Parameters
    ----------
    plan : str, optional
        ["A", "B"], Default = "A"
        * Plan A is for a fully populated mirror (798 segments)
        * Plan B has the inner 5 rings missing (588 segments) and a further
        5 random segments missing (583 segments)
    wavelength : float, list, array, optional
        [um] Default = 2.2um. The wavelength(s) for which a PSF should be made
    mode : str, optional
        ["wide", "zoom"] Default = "wide". Sets the pixel size for each
        of the MICADO imaging modes - {"wide" : 4mas, "zoom" : 1.5mas}
    size : int, optional
        [pixels] Default = 1024
    segments : list, optional
        Default = None. A list of which segments to use for generating the E-ELT
        mirror. See ``get_eelt_segments()``
    filename : str, optional
        Default = None. If filename is not None, the resulting FITS object will
        be saved to disk
    use_pupil_mask : str, optional
        Default = True.

    Other Parameters
    ----------------
    Values to pass to the POPPY functions

    flattoflat : float
        [m] Default : 1.256
    gap  : float
        [m] Default  : 0.004
    secondary_radius  : float
        [m] Default : 5
    n_supports : int
        Default : 6
    support_width : float
        [m] Default : 0.2
    support_angle_offset  : float
        [deg] Default : 0
    n_missing : int
        Default : None. Number of segments missing
    pupil_inner_radius  : float
        [m] Default : None  # Plan A: 5.6m, Plan B: 11.5m
    pupil_outer_radius  : float
        [m] Default : 19


    Returns
    -------
    ``astropy.HDUList`` : an astropy FITS object with the PSF in the data
    extensions


    See also
    --------
    :func:`.get_eelt_segments`

    """

    try:
        import poppy
    except:
        raise ImportError("""Poppy is not installed -
        See https://pythonhosted.org/poppy""")

    params = {"flattoflat"           : 1.256,
              "gap"                  : 0.004,
              "secondary_radius"     : 5,
              "n_supports"           : 6,
              "support_width"        : 0.2,
              "support_angle_offset" : 0,
              "n_missing"            : None,
              "oversample"           : 2,
              "pupil_inner_radius"    : None,
              "pupil_outer_radius"    : 19}
    # Careful - when calling the detector, the rusulting PSF is 2x oversampled.
    # Hence we double the pixelscale at osys.add_detector

    params.update(**kwargs)

    params["pixelscale"] = 0.004 if mode.lower() == "wide" else 0.0015
    if params["pupil_inner_radius"] is None:
        params["pupil_inner_radius"] = 11.5 if plan.lower() == "b" else 5.6



    if segments is None:
        segments = get_eelt_segments(plan=plan, missing=params["n_missing"])

    ap = poppy.MultiHexagonAperture(flattoflat = params["flattoflat"],
                                    gap        = params["gap"],
                                    segmentlist= segments)
    sec = poppy.SecondaryObscuration(secondary_radius    = params["secondary_radius"],
                                     n_supports          = params["n_supports"],
                                     support_width       = params["support_width"],
                                     support_angle_offset= params["support_angle_offset"])

    opticslist = [ap, sec]

    if use_pupil_mask:
        cold_in = poppy.SecondaryObscuration(secondary_radius=params["pupil_inner_radius"], n_supports=0)
        cold_out = poppy.CircularAperture(radius=params["pupil_outer_radius"])
        opticslist += [cold_in, cold_out]

    eelt = poppy.CompoundAnalyticOptic(opticslist=opticslist,
                                       name='E-ELT Plan '+plan)

    osys = poppy.OpticalSystem(oversample=params["oversample"], pupil_diameter=50)
    osys.add_pupil(eelt)
    osys.add_detector(pixelscale = params["pixelscale"],
                      fov_arcsec = params["pixelscale"] * size)

    if np.any(wavelength) < 0.1:
        logging.warning("One or more wavelengths is/are very short")
        print(wavelength)

    wavelength *= 1E-6
    hdu_list = osys.calc_psf(wavelength)
    hdu_list[0].data = spi.zoom(hdu_list[0].data, 1./params["oversample"])
    hdu_list[0].data /= np.sum(hdu_list[0].data)

    for hdu in hdu_list:
        hdu.data = hdu.data.astype(np.float32)
        hdu.header["PIXELSCL"] = params["pixelscale"]
        hdu.header["PSF_TYPE"] = "POPPY"
        hdu.header["CDELT1"]   = params["pixelscale"]
        hdu.header["CDELT2"]   = params["pixelscale"]
        hdu.header["CUNIT1"]   = "arcsec"
        hdu.header["CUNIT2"]   = "arcsec"

    if filename is None:
        return hdu_list
    else:
        print("Writing to", filename)
        hdu_list.writeto(filename, clobber=True)



def get_eelt_segments(plan="A", missing=None, return_missing_segs=False,
                      inner_diam=10.6, outer_diam=39.):
    """
    Generate a list of segments for POPPY for the E-ELT

    Parameters
    ----------
    plan : str, optional
        ["A", "B"], Default = "A"
        * Plan A is for a fully populated mirror (798 segments)
        * Plan B has the inner 5 rings missing (588 segments) and a further
        5 random segments missing (583 segments)
    missing : int, list, optional
        Default = None. If an integer is passed, this many random segments are
        removed. If ``missing`` is a list, entries refer to specific segment IDs
    return_missing_segs : bool, optional
        Defualt is False. Returns the missing segment numbers
    inner_diam : float, optional
        [m] Default = 10.6. Diameter which produces ESO's mirror configuration
    outer_diam : float, optional
        [m] Default = 39.0. Diameter which produces ESO's mirror configuration


    Returns
    -------
    segs : list
        A list of segment IDs for the mirror segments.
    missing : list, conditional
        Only returned if ``return_missing_segs == True``. A list of segment IDS
        for the segments which are missing.

    """

    try:
        import poppy
    except:
        raise ImportError("Poppy is not installed - google 'JWST POPPY'")

    if plan.lower() == "b":
        #inner_diam = 21.9
        first_seg = 270
        if missing is None: missing = 5
    else:
        first_seg = 60
        if missing is None: missing = 0


    ap = poppy.MultiHexagonAperture(flattoflat=1.256, gap=0.004,
                                                    segmentlist=np.arange(2000))
    rad = [np.sqrt(np.sum(np.array(ap._hex_center(j))**2))
                                                        for j in ap.segmentlist]

    mask = (np.array(rad) < 0.5*outer_diam) * (ap.segmentlist > first_seg)
    segs = np.array(ap.segmentlist)[mask]

    if isinstance(missing, int):
        if missing > 0:
            i = np.random.randint(0, len(segs)-missing, size=missing)
            missing = segs[i]
        else:
            missing = []

    for i in missing:
        if i in segs:
            x = np.where(segs == i)[0][0]
            segs = segs.tolist()
            segs.pop(x)
            segs = np.array(segs)

    if return_missing_segs:
        return segs, missing
    else:
        return segs
