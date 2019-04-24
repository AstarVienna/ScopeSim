# pylint: disable=too-many-lines
"""
The module that contains the functionality to create Source py_objects

Module Summary
--------------
The Source is essentially a list of spectra and a list of positions. The
list of positions contains a reference to the relevant spectra. The advantage
here is that if there are repeated spectra in a data cube, we can reduce the
amount of calculations needed. Furthermore, if the input is originally a list
of stars, etc., where the position of a star is not always an integer multiple
of the plate scale, we can keep the information until the PSFs are needed.

Classes
-------
Source

Functions
---------
Functions to create ``Source`` py_objects
::

    empty_sky()
    star(mag, filter_name="Ks", spec_type="A0V", x=0, y=0)
    stars(mags, filter_name="Ks", spec_types=["A0V"], x=[0], y=[0])
    star_grid(n, mag_min, mag_max, filter_name="Ks", separation=1, area=1,
              spec_type="A0V")
    source_from_image(images, lam, spectra, pix_res, oversample=1,
                      units="ph/s/m2", flux_threshold=0,
                      center_pixel_offset=(0, 0))
    source_1E4_Msun_cluster(distance=50000, half_light_radius=1)


Functions for manipulating spectra for a ``Source`` object
::

    scale_spectrum(lam, spec, mag, filter_name="Ks", return_ec=False)
    scale_spectrum_sb(lam, spec, mag_per_arcsec, pix_res=0.004,
                      filter_name="Ks", return_ec=False)
    flat_spectrum(mag, filter_name="Ks", return_ec=False)
    flat_spectrum_sb(mag_per_arcsec, filter_name="Ks", pix_res=0.004,
                     return_ec=False)


Functions regarding photon flux and magnitudes
::

    zero_magnitude_photon_flux(filter_name)
    _get_stellar_properties(spec_type, cat=None, verbose=False)
    _get_stellar_mass(spec_type)
    _get_stellar_Mv(spec_type)
    _get_pickles_curve(spec_type, cat=None, verbose=False)


Helper functions
::

    value_at_lambda(lam_i, lam, val, return_index=False)
    SED(spec_type, filter_name="V", magnitude=0.)

"""
###############################################################################
# source
#
# DESCRIPTION

#
# The source contains two arrays:
#  - PositionArray:
#  - SpectrumArray


# Flow of events
# - Generate the lists of spectra and positions
# - Apply the transmission curves [SpectrumArray]
# - shrink the 1D spectra to the resolution of the psf layers [SpectrumArray]
# - apply any 2D spatial effects [PositionArray]
# for i in len(slices)
#   - generate a working slice [PositionArray, SpectrumArray, WorkingSlice]
#   - apply the PSF for the appropriate wavelength [WorkingSlice]
#   - apply any wavelength dependent spatials [WorkingSlice]
#   - apply Poisson noise to the photons in the slice [WorkingSlice]
#   - add the WorkingSlice to the FPA [WorkingSlice, FPArray]


# TODO implement conversions to Source object from:
# ascii
#    x, y, mag, [temp]
#    x, y, type
# images
#    JHK
#    cube


import os

from copy import deepcopy

import numpy as np
from scipy.ndimage import sum as ndisum
import scipy.ndimage.interpolation as spi
from scipy.signal import fftconvolve

from astropy.io import fits
from astropy.convolution import convolve
import astropy.units as u

from .OLD_source_utils import scale_spectrum, spectrum_sum_over_range
from ..OLD_spectral import TransmissionCurve, EmissionCurve,\
    UnityCurve, BlackbodyCurve
from scopesim.optics import OLD_psf as sim_psf
from .. import utils


# add_uniform_background() moved to detector
# get_slice_photons() renamed to photons_in_range() and moved to Source
# _apply_transmission_curve() moved to Source
# apply_optical_train() moved to Source


class Source(object):
    """
    Create a source object from a file or from arrays

    Source class generates the arrays needed for source. It takes various
    inputs and converts them to an array of positions and references to spectra.
    It also converts spectra to photons/s/voxel. The default units for input
    data is ph/s/m2/bin.

    The internal variables are related like so:
    ::
        f(x[i], y[i]) = spectra[ref[i]] * weight[i]


    Parameters
    ----------
    filename : str
        FITS file that contains either a previously saved ``Source`` object or a
        data cube with dimensions x, y, lambda. A ``Source`` object is
        identified by the header keyword SIMCADO with value SOURCE.

    or

    lam : np.array
        [um] Wavelength bins of length (m)
    spectra : np.array
        [ph/s/m2/bin] A (n, m) array with n spectra, each with m spectral bins
    x, y : np.array
        [arcsec] coordinates of where the emitting files are relative to the
        centre of the field of view
    ref : np.array
        the index for .spectra which connects a position (x, y) to a spectrum
        f(x[i], y[i]) = spectra[ref[i]] * weight[i]
    weight : np.array
        A weighting to scale the relevant spectrum for each position

    Keyword arguments
    -----------------
    units : str
        The units of the spectra. Default is ph/s/m2/bin
    pix_unit : str
        Default is "arcsec". Acceptable are "arcsec", "arcmin", "deg", "pixel"
    exptime : float
        If the input spectrum is not normalised to 1 sec
    area : float
        The telescope area used to generate the source object
    pix_res : float
        [arcsec] The pixel resolution of the detector. Useful for surface
        brightness calculations
    bg_spectrum : EmissionCurve
        If there is a surface brightness term to add, add it here

    """

    def __init__(self, filename=None,
                 lam=None, spectra=None, x=None, y=None, ref=None, weight=None,
                 **kwargs):

        self.params = {"units"   : "ph/s",
                       "pix_unit": "arcsec",
                       "exptime" : 1,
                       "area"    : 1,
                       "pix_res" : 0.004,
                       "bg_spectrum" : None}
        self.params.update(kwargs)

        if isinstance(x, (tuple, list)):
            x = np.array(x)
        if isinstance(y, (tuple, list)):
            y = np.array(y)

        if x is not None:
            x = x.astype(np.float32)
        if y is not None:
            y = y.astype(np.float32)

        if "pix" in self.params["pix_unit"]:
            x *= self.params["pix_res"]
            y *= self.params["pix_res"]
        elif "arcmin" in self.params["pix_unit"]:
            x *= 60.
            y *= 60.
        elif "deg" in self.params["pix_unit"]:
            x *= 3600.
            y *= 3600.

        self.info = {}
        self.info['description'] = "List of spectra and their positions"

        self.units = u.Unit(self.params["units"])
        self.exptime = self.params["exptime"]
        self.pix_res = self.params["pix_res"]

        self.x = None
        self.y = None  # set later

        # A file can contain a previously saved Source object; in this case the
        # header keyword
        # "SIMCADO" is set to "SOURCE". If this is not the case, we assume that
        # the file
        # contains a data cube with dimensions x, y, lambda.
        # If no filename is given, we build the Source from the arrays.
        if filename is not None:
            hdr = fits.getheader(filename)
            if "SIMCADO" in hdr.keys() and hdr["SIMCADO"] == "SOURCE":
                self.read(filename)
            else:
                self._from_cube(filename)
        elif not any(elem is None for elem in (lam, spectra, x, y, ref)):
            self._from_arrays(lam, spectra, x, y, ref, weight)
        else:
            raise ValueError("Trouble with inputs. Could not create Source")

        self.ref = np.array(self.ref, dtype=int)
        self.x_orig = deepcopy(self.x)
        self.y_orig = deepcopy(self.y)
        self.spectra_orig = deepcopy(self.spectra)

        self.bg_spectrum = None

    @classmethod
    def load(cls, filename):
        """Load :class:'.Source' object from filename"""
        import pickle
        with open(filename, 'rb') as fp1:
            src = pickle.load(fp1)
        return src

    def dump(self, filename):
        """Save to filename as a pickle"""
        import pickle
        with open(filename, 'wb') as fp1:
            pickle.dump(self, fp1)

    def apply_optical_train(self, opt_train, detector, chips="all",
                            sub_pixel=False, **kwargs):
        """
        Apply all effects along the optical path to the source photons

        Parameters
        ----------
        opt_train : scopesim.OpticalTrain
            the object containing all information on what is to happen to the
            photons as they travel from the source to the detector
        detector : scopesim.Detector
            the object representing the detector
        chips : int, str, list, optional
            The IDs of the chips to be readout. "all" is also acceptable
        sub_pixel : bool, optional
            if sub-pixel accuracy is needed, each source is shifted individually.
            Default is False

        Other Parameters
        ----------------
        INST_DEROT_PERFORMANCE : float
            [0 .. 100] Percentage of the sky rotation that the derotator removes
        SCOPE_JITTER_FWHM : float
            [arcsec] The FWMH of the gaussian blur caused by jitter
        SCOPE_DRIFT_DISTANCE : float
            [arcsec] How far from the centre of the field of view has the
            telescope drifted during a DIT

        Notes
        -----
        Output array is in units of [ph/s/pixel] where the pixel is internal
        oversampled pixels - not the pixel size of the detector chips

        """
        params = {"verbose"                : opt_train.cmds.verbose,
                  "INST_DEROT_PERFORMANCE" : opt_train.cmds["INST_DEROT_PERFORMANCE"],
                  "SCOPE_JITTER_FWHM"      : opt_train.cmds["SCOPE_JITTER_FWHM"],
                  "SCOPE_DRIFT_DISTANCE"   : opt_train.cmds["SCOPE_DRIFT_DISTANCE"],
                  "sub_pixel"              : sub_pixel}
        params.update(self.params)
        params.update(kwargs)

        self.pix_res = opt_train.pix_res

        # 1. Apply the master transmission curve to all the spectra
        #
        # 1.5 Create a canvas onto which we splat the PSFed files
        #
        # 2. For each layer between cmds.lam_bin_edges[i, i+1]
        #   - Apply the x,y shift for the ADC
        #       - Apply any other shifts
        #   - Apply the PSF for the layer
        #       - Sum up the photons in this section of the spectra
        #   - Add the layer to the final image array
        #
        # 3. Apply wave-indep psfs
        #   - field rotation
        #   - telescope shake
        #   - tracking error
        #
        # 3.5 Up until now everything is ph/s/m2/bin
        #     Apply the collecting area of the telescope
        #
        # 4. Add the average number of atmo-bg and mirror-bb photons
        # 5. Apply the instrumental distortion

        if chips is None or str(chips).lower() == "all":
            chips = np.arange(len(detector.chips))

        if not hasattr(chips, "__len__"):
            chips = [chips]

        # 1.
        self._apply_transmission_curve(opt_train.tc_source)

        for chip_i in chips:
            print("Generating image for chip", detector.chips[chip_i].id)

            # 1.5
            image = None

            # 2.
            for i in range(len(opt_train.lam_bin_edges[:-1])):

                if params["verbose"]:
                    print("Wavelength slice [um]:",
                          opt_train.lam_bin_centers[i])

                # apply the adc shifts
                self._x = self.x + opt_train.adc_shifts[0][i]
                self._y = self.y + opt_train.adc_shifts[1][i]

                # include any other shifts here

                # apply the psf (get_slice_photons is called within)
                lam_min, lam_max = opt_train.lam_bin_edges[i:i+2]
                psf_i = utils.nearest(opt_train.psf.lam_bin_centers,
                                      opt_train.lam_bin_centers[i])
                psf = opt_train.psf[psf_i]

                oversample = opt_train.cmds["SIM_OVERSAMPLING"]
                sub_pixel = params["sub_pixel"]
                verbose = params["verbose"]

                # image is in units of ph/s/pixel/m2
                imgslice = self.image_in_range(psf, lam_min, lam_max,
                                               detector.chips[chip_i],
                                               pix_res=opt_train.pix_res,
                                               oversample=oversample,
                                               sub_pixel=sub_pixel,
                                               verbose=verbose)
                if image is None:
                    image = imgslice
                else:
                    image += imgslice

            # 3. Apply wavelength-independent spatial effects
            # !!!!!!!!!!!!!! All of these need to be combined into a single
            # function that traces out the path taken by the telescope,
            # rather than having the arcs from the derotator() function
            # being stretched by the tracking() function and then the whole
            # thing blurred by wind_jitter()
            if params["INST_DEROT_PERFORMANCE"] < 100:
                image = opt_train.apply_derotator(image)
            if params["SCOPE_DRIFT_DISTANCE"] > 0.33 * self.pix_res:
                image = opt_train.apply_tracking(image)
            if params["SCOPE_JITTER_FWHM"] > 0.33 * self.pix_res:
                image = opt_train.apply_wind_jitter(image)

            # 3.5 Scale by telescope area
            image *= opt_train.cmds.area

            # 4. Add backgrounds
            image += (opt_train.n_ph_atmo + opt_train.n_ph_mirror +
                      opt_train.n_ph_ao)

            # TODO: protected members should not be set by another class (OC)
            #       These could be added to info dictionary, if they're only
            #       informational.
            detector._n_ph_atmo = opt_train.n_ph_atmo
            detector._n_ph_mirror = opt_train.n_ph_mirror
            detector._n_ph_ao = opt_train.n_ph_ao

            # 5. Project onto chip
            self.project_onto_chip(image, detector.chips[chip_i])

        ######################################
        # CAUTION WITH THE PSF NORMALISATION #
        ######################################

    def project_onto_chip(self, image, chip):
        """
        Re-project the photons onto the same grid as the detectors use

        Parameters
        ----------
        image : np.ndarray
            the image to be re-projected
        chip : detector.Chip
            the chip object where the image will land
        """
        # This is just a change of pixel scale
        chip.reset()
        scale_factor = self.pix_res / chip.pix_res
        chip_arr = spi.zoom(image, scale_factor, order=1)
        chip_arr *= np.sum(image) / np.sum(chip_arr)
        chip.add_signal(chip_arr)

    def image_in_range(self, psf, lam_min, lam_max, chip, **kwargs):
        """
        Find the files that fall in the chip area and generate an image for
        the wavelength range [lam_min, lam_max)

        Output is in [ph/s/pixel]

        Parameters
        ----------
        psf : psf.PSF object
            The PSF that the files will be convolved with
        lam_min, lam_max : float
            [um] the wavelength range relevant for the psf
        chip : str, detector.Chip
            - detector.Chip : the chip that will be seeing this image.
            - str : ["tiny", "small", "center"] -> [128, 1024, 4096] pixel chips


        Optional parameters (**kwargs)
        ------------------------------
        sub_pixel : bool
            if sub-pixel accuracy is needed, each source is shifted individually
            Default is False
        pix_res : float
            [arcsec] the field of view of each pixel. Default is 0.004 arcsec
        oversample : int
            the psf images will be oversampled to better conserve flux.
            Default is 1 (i.e. not oversampled)
        verbose : bool
            Default that of the OpticalTrain object

        Returns
        -------
        slice_array : np.ndarray
            the image of the source convolved with the PSF for the given range

        """

        params = {"pix_res"     : 0.004,
                  "sub_pixel"   : False,
                  "oversample"  : 1,
                  "verbose"     : False}

        params.update(kwargs)

        #  no PSF given: use a delta kernel
        if isinstance(psf, type(None)):
            psf = np.zeros((7, 7))
            psf[3, 3] = 1

        # psf cube given: extract layer for central wavelength
        if isinstance(psf, (sim_psf.PSFCube, sim_psf.UserPSFCube)):
            lam_cen = (lam_max + lam_min) / 2.
            ##############################################################
            # Bad coding - psf changes from PSFCube object to PSF object #
            ##############################################################
            # ..todo :: fix this
            psf = psf.nearest(lam_cen)

        # psf given as array: convert to PSF object
        if isinstance(psf, np.ndarray):
            arr = deepcopy(psf)
            pix_res = params["pix_res"] / params["oversample"]
            size = psf.shape[0]
            psf = sim_psf.PSF(size, pix_res)
            psf.set_array(arr)

        # .. TODO: There is no provision for chip rotation wrt (x,y) system (OC)
        # Create Chip object if chip described by a string
        if isinstance(chip, str):
            from ..detector import Chip
            if chip.lower() == "small":
                chip = Chip(0, 0, 1024, 1024, 0.004)
            elif "cent" in chip.lower():
                chip = Chip(0, 0, 4096, 4096, 0.004)
            elif "tiny" in chip.lower():
                chip = Chip(0, 0, 128, 128, 0.004)
            else:
                raise ValueError("Unknown chip identification")

        # Check whether _x has been created - _x contains the adc corrections
        if not hasattr(self, "_x"):
            self._x = np.copy(self.x)
            self._y = np.copy(self.y)

        # Determine x- and y- range covered by chip
        # TODO: Use chip.wcs to convert (x, y) into pixel coordinates,
        #       then simply cut at the pixel edges. Alternatively,
        #       project chip edges to the sky.
        if chip is not None:
            mask = (self._x > chip.x_min) * (self._x < chip.x_max) * \
                   (self._y > chip.y_min) * (self._y < chip.y_max)
            params["pix_res"] = chip.pix_res / params["oversample"]
            x_min, x_max = chip.x_min, chip.x_max,
            y_min, y_max = chip.y_min, chip.y_max
            x_cen, y_cen = chip.x_cen, chip.y_cen

            naxis1, naxis2 = chip.naxis1, chip.naxis2

        else:
            # no chip given: use area covered by object arrays
            mask = np.array([True] * len(self._x))
            params["pix_res"] /= params["oversample"]
            x_min, x_max = np.min(self._x), np.max(self._x)
            y_min, y_max = np.min(self._y), np.max(self._y),
            x_cen, y_cen = (x_max + x_min) / 2, (y_max + y_min) / 2

            # the conversion to int was causing problems because some
            # values were coming out at 4095.9999, so the array was (4095, 4096)
            # hence the 1E-3 on the end
            naxis1 = int((x_max - x_min) / params["pix_res"] + 1E-3)
            naxis2 = int((y_max - y_min) / params["pix_res"] + 1E-3)

        slice_array = np.zeros((naxis1, naxis2), dtype=np.float32)
        slice_photons = self.photons_in_range(lam_min, lam_max)

        # convert point source coordinates to pixels
        x_pix = (self._x - x_cen) / params["pix_res"]
        y_pix = (self._y - y_cen) / params["pix_res"]

        self.x_pix = x_pix + chip.naxis1 // 2
        self.y_pix = y_pix + chip.naxis2 // 2

        # if sub-pixel accuracy is needed, be prepared to wait. For this we
        # need to go through every source spectrum in turn, shift the psf by
        # the decimal amount given by pos - int(pos), then place a
        # certain slice of the psf on the output array.
        ax, ay = np.array(slice_array.shape) // 2
        bx, by = np.array(psf.array.shape)   // 2
        mx, my = np.array(psf.array.shape) % 2

        if params["verbose"]:
            print("Chip ID:", chip.id,
                  "- Creating layer between [um]:", lam_min, lam_max)

        psf_array = np.copy(psf.array)

        if params["sub_pixel"] is True:
            # for each point source in the list, add a psf to the slice_array
            # x_int, y_int = np.floor(x_pix), np.floor(y_pix)
            # dx, dy = src.x - x_int, src.y - y_int

            if bx == ax and by == ay:
                pass
            elif bx > ax and by > ay:
                # psf_array larger than slice_array: cut down
                psf_array = psf_array[(bx - ax):(bx + ax), (by - ay):(by + ay)]
            elif bx < ax and by < ay:
                # psf_array smaller than slice_array: pad with zeros
                pad_x, pad_y = ax - bx, ay - by
                psf_array = np.pad(psf_array,
                                   ((pad_x, pad_x-mx),
                                    (pad_y, pad_y-my)),
                                   mode="constant")
            else:
                print("PSF", psf.array.shape, "Chip", slice_array.shape)
                raise ValueError("PSF and Detector chip sizes are odd:")

            for i in range(len(x_pix)):
                psf_tmp = np.copy(psf_array)
                print(x_pix[i], y_pix[i])
                psf_tmp = spi.shift(psf_tmp, (x_pix[i], y_pix[i]), order=1)
                slice_array += psf_tmp * slice_photons[i]

        elif params["sub_pixel"] == "raw":
            x_int, y_int = np.floor(x_pix), np.floor(y_pix)
            i = (ax + x_int[mask]).astype(int)
            j = (ay + y_int[mask]).astype(int)
            slice_array[i, j] = slice_photons[mask]

        else:
            # If astrometric precision is not that important and everything
            # has been oversampled, use this section.
            #  - ax, ay are the pixel coordinates of the image centre
            # use np.floor instead of int-ing
            x_int, y_int = np.floor(x_pix), np.floor(y_pix)
            i = (ax + x_int[mask]).astype(int)
            j = (ay + y_int[mask]).astype(int)

            # The following is faster than a loop
            ij = i * naxis1 + j   # naxis1 or naxis2?
            iju = np.unique(ij)
            slice_array.flat[iju] += ndisum(slice_photons[mask].flat,
                                            ij, iju)

            try:
                # slice_array = convolve_fft(slice_array, psf.array,
                #                            allow_huge=True)
                # make the move to scipy
                slice_array = fftconvolve(slice_array, psf.array, mode="same")
            except ValueError:
                slice_array = convolve(slice_array, psf.array)

        return slice_array

    def photons_in_range(self, lam_min=None, lam_max=None):
        """

        Number of photons between lam_min and lam_max in units of [ph/s/m2]

        Calculate how many photons for each source exist in the wavelength range
        defined by lam_min and lam_max.

        Parameters
        ----------
        lam_min, lam_max : float, optional
            [um] integrate photons between these two limits. If both are ``None``,
            limits are set at lam[0], lam[-1] for the source's wavelength range

        Returns
        -------
        slice_photons : float
            [ph/s/m2] The number of photons in the wavelength range

        """
        spec_photons = spectrum_sum_over_range(self.lam, self.spectra, lam_min, lam_max)

        slice_photons = spec_photons[self.ref] * self.weight
        return slice_photons

    def scale_spectrum(self, idx=0, mag=20, filter_name="Ks"):
        """
        Scale a certain spectrum to a certain magnitude

        See :func:`scopesim.source.scale_spectrum` for examples

        Parameters
        ----------
        idx : int
            The index of the spectrum to be scaled: <Source>.spectra[idx]
            Default is <Source>.spectra[0]
        mag : float
            [mag] new magnitude of spectrum
        filter_name : str, TransmissionCurve
           Any filter name from SimCADO or a
           :class:`~.scopesim.spectral.TransmissionCurve` object
           (see :func:`~.scopesim.optics.get_filter_set`)
        """

        self.lam, self.spectra[idx] = scale_spectrum(lam=self.lam,
                                                     spec=self.spectra[idx],
                                                     mag=mag,
                                                     filter_name=filter_name,
                                                     return_ec=False)


    def scale_with_distance(self, distance_factor):
        """
        Scale the source for a new distance

        Scales the positions and brightnesses of the :class:`.Source` object
        according to the ratio of the new and old distances

        i.e. distance_factor = new_distance / current_distance

        .. warning::
            This does not yet take into account redshift

        .. todo::
            Implement redshift

        Parameters
        ----------
        distance_factor : float
            The ratio of the new distance to the current distance
            i.e. distance_factor = new_distance / current_distance

        Examples
        --------
        ::

            >>> from scopesim.source.templates import cluster
            >>>
            >>> curr_dist = 50000  # pc, i.e. LMC
            >>> new_dist = 770000  # pc, i.e. M31
            >>> src = cluster(distance=curr_dist)
            >>> src.scale_with_distance( new_dist/curr_dist )

        """
        self.x /= distance_factor
        self.y /= distance_factor
        self.weight /= distance_factor**2


    def add_background_surface_brightness(self):
        """
        Add an EmissionCurve for the background surface brightness of the object
        """
        pass


    def rotate(self, angle, unit="degree", use_orig_xy=False):
        """
        Rotates the ``x`` and ``y`` coordinates by ``angle`` [degrees]

        Parameters
        ----------
        angle : float
            Default is in degrees, this can set with ``unit``
        unit : str, astropy.Unit
            Either a string with the unit name, or an
            ``astropy.unit.Unit`` object
        use_orig_xy : bool
            If the rotation should be based on the original coordinates or the
            current coordinates (e.g. if rotation has already been applied)

        """
        ang = (angle * u.Unit(unit)).to(u.rad)

        if use_orig_xy:
            xold, yold = self.x_orig, self.y_orig
        else:
            xold, yold = self.x, self.y

        self.x = xold * np.cos(ang) - yold * np.sin(ang)
        self.y = xold * np.sin(ang) + yold * np.cos(ang)


    def shift(self, dx=0, dy=0, use_orig_xy=False):
        """
        Shifts the coordinates of the source by (dx, dy) in [arcsec]

        Parameters
        ----------
        dx, dy : float, array
            [arcsec] The offsets for each coordinate in the arrays ``x``, ``y``.
            - If dx, dy are floats, the same offset is applied to all coordinates
            - If dx, dy are arrays, they must be the same length as ``x``, ``y``
        use_orig_xy : bool
            If the shift should be based on the original coordinates or the
            current coordinates (e.g. if shift has already been applied)

        """
        self.dx = dx
        self.dy = dy

        if use_orig_xy:
            self.x = self.x_orig + dx
            self.y = self.y_orig + dy
        else:
            self.x += dx
            self.y += dy


    def on_grid(self, pix_res=0.004):
        """
        Return an image with the positions of all files.

        The pixel values correspond to the number of emitting py_objects in that
        pixel

        Parameters
        ----------
        pix_res : float
            [arcsec] The grid spacing

        Returns
        -------
        im : 2D array
            A numpy array containing an image of where the files are

        """

        xmin = np.min(self.x)
        ymin = np.min(self.y)
        x_i = ((self.x - xmin) / pix_res).astype(int)
        y_i = ((self.y - ymin) / pix_res).astype(int)
        img = np.zeros((np.max(x_i)+2, np.max(y_i)+2))
        img[x_i, y_i] += 1

        return img


    def read(self, filename):
        """
        Read in a previously saved :class:`.Source` FITS file

        Parameters
        ----------
        filename : str
            Path to the file

        """

        ipt = fits.open(filename)
        dat0 = ipt[0].data
        hdr0 = ipt[0].header
        dat1 = ipt[1].data
        hdr1 = ipt[1].header
        ipt.close()

        self.x = dat0[0, :]
        self.y = dat0[1, :]
        self.ref = dat0[2, :]
        self.weight = dat0[3, :]

        lam_min, lam_max = hdr1["LAM_MIN"], hdr1["LAM_MAX"]
        self.lam_res = hdr1["LAM_RES"]
        self.lam = np.linspace(lam_min, lam_max, hdr1["NAXIS1"])
        self.spectra = dat1

        if "BUNIT" in hdr0.keys():
            self.params["units"] = u.Unit(hdr0["BUNIT"])
        if "EXPTIME" in hdr0.keys():
            self.params["exptime"] = hdr0["EXPTIME"]
        if "AREA"   in hdr0.keys():
            self.params["area"] = hdr0["AREA"]
        if "CDELT1" in hdr0.keys():
            self.params["pix_res"] = hdr0["CDELT1"]
        if "CUNIT1" in hdr0.keys():
            self.params["pix_unit"] = u.Unit(hdr0["CUNIT1"])
        self.lam_res = hdr1["LAM_RES"]

        self._convert_to_photons()


    def write(self, filename):
        """
        Write the current Source object out to a FITS file

        Parameters
        ----------
        filename : str
            where to save the FITS file

        Notes
        -----
        Just a place holder so that I know what's going on with the input table
        * The first extension [0] contains an "image" of size 4 x N where N is the
        number of files. The 4 columns are x, y, ref, weight.
        * The second extension [1] contains an "image" with the spectra of all
        files. The image is M x len(spectrum), where M is the number of unique
        spectra in the source list. M = max(ref) - 1
        """

        # hdr = fits.getheader("../../../PreSim/Input_cubes/GC2.fits")
        # ipt = fits.getdata("../../../PreSim/Input_cubes/GC2.fits")
        # flux_map = np.sum(ipt, axis=0).astype(dtype=np.float32)
        # x,y = np.where(flux_map != 0)
        # ref = np.arange(len(x))
        # weight = np.ones(len(x))
        # spectra = np.swapaxes(ipt[:,x,y], 0, 1)
        # lam = np.linspace(0.2,2.5,231)

        xyHDU = fits.PrimaryHDU(np.array((self.x, self.y, self.ref, self.weight)))
        xyHDU.header["X_COL"] = "1"
        xyHDU.header["Y_COL"] = "2"
        xyHDU.header["REF_COL"] = "3"
        xyHDU.header["W_COL"] = "4"

        xyHDU.header["BUNIT"] = self.units.to_string()
        xyHDU.header["EXPTIME"] = self.params["exptime"]
        xyHDU.header["AREA"] = self.params["area"]
        xyHDU.header["CDELT1"] = self.params["pix_res"]
        xyHDU.header["CDELT2"] = self.params["pix_res"]
        xyHDU.header["CUNIT1"] = self.params["pix_unit"]
        xyHDU.header["CUNIT2"] = self.params["pix_unit"]

        xyHDU.header["SIMCADO"] = "SOURCE"

        specHDU = fits.ImageHDU(self.spectra)
        specHDU.header["CRVAL1"] = self.lam[0]
        specHDU.header["CRPIX1"] = 0
        specHDU.header["CDELT1"] = (self.lam_res, "[um] Spectral resolution")
        specHDU.header["LAM_MIN"] = (self.lam[0], "[um] Minimum wavelength")
        specHDU.header["LAM_MAX"] = (self.lam[-1], "[um] Maximum wavelength")
        specHDU.header["LAM_RES"] = (self.lam_res, "[um] Spectral resolution")

        hdu = fits.HDUList([xyHDU, specHDU])
        hdu.writeto(filename, overwrite=True, checksum=True)


    @property
    def info_keys(self):
        """Return keys of the `info` dict"""
        return self.info.keys()


    def _apply_transmission_curve(self, transmission_curve):
        """
        Apply the values from a TransmissionCurve object to self.spectra

        Parameters
        ----------
        transmission_curve : TransmissionCurve
            The TransmissionCurve to be applied

        See Also
        --------
        :class:`scopesim.spectral.TransmissionCurve`

        """
        tc = deepcopy(transmission_curve)
        tc.resample(self.lam, use_default_lam=False)
        self.spectra = self.spectra_orig * tc.val


    def _convert_to_photons(self):
        """
        Convert the spectra to photons/(s m2)

        If [arcsec] are in the units, we want to find the photons per pixel.
        If [um] are in the units, we want to find the photons per wavelength bin.

        .. todo::
            Come back and put in other energy units like Jy, mag, ergs

        """

        self.units = u.Unit(self.params["units"])
        bases = self.units.bases

        factor = u.Quantity(1.)
        if u.s not in bases:
            factor /= (self.params["exptime"] * u.s)
        if u.m not in bases:
            factor /= (1. * u.m**2)
        if u.micron in bases:
            factor *= (self.lam_res * u.um)
        if u.arcsec in bases:
            factor *= (self.params["pix_res"] * u.arcsec)**2

        self.units = self.units * factor.unit
        self.spectra *= factor.value


    def _from_cube(self, filename):
        # Should this be a class method?
        """
        Make a Source object from a cube in memory or a FITS cube on disk

        Parameters
        ----------
        filename : str
            Path to the FITS cube

        """

        if isinstance(filename, str) and os.path.exists(filename):
            hdr = fits.getheader(filename)
            cube = fits.getdata(filename)
        else:
            raise ValueError(filename + " doesn't exist")

        lam_res = hdr["CDELT3"]
        lam_min = hdr["CRVAL3"] - hdr["CRPIX3"] * lam_res
        lam_max = lam_min + hdr["NAXIS3"] * lam_res

        flux_map = np.sum(cube, axis=0).astype(dtype=np.float32)
        x, y = np.where(flux_map != 0)

        self.lam = np.linspace(lam_min, lam_max, hdr["NAXIS3"])
        self.spectra = np.swapaxes(cube[:, x, y], 0, 1)
        self.x = x
        self.y = y
        self.ref = np.arange(len(x))
        self.weight = np.ones(len(x))

        if "BUNIT" in hdr.keys():
            self.params["units"] = u.Unit(hdr["BUNIT"])
        if "EXPTIME" in hdr.keys():
            self.params["exptime"] = hdr["EXPTIME"]
        if "AREA"   in hdr.keys():
            self.params["area"] = hdr["AREA"]
        if "CDELT1" in hdr.keys():
            self.params["pix_res"] = hdr["CDELT1"]
        if "CUNIT1" in hdr.keys():
            self.params["pix_unit"] = hdr["CUNIT1"]
        self.lam_res = lam_res

        self._convert_to_photons()


    def _from_arrays(self, lam, spectra, x, y, ref, weight=None):
        # Should this be a class method?
        """
        Make a Source object from a series of lists

        Parameters
        ----------
        lam : np.ndarray
            Dimensions (1, m) with m spectral bins
        spectra : np.ndarray
            Dimensions (n, m) for n SEDs, each with m spectral bins
        x, y : np.ndarray
            [arcsec] each (1, n) for the coordinates of n emitting py_objects
        ref : np.ndarray
            Dimensions (1, n) for referencing each coordinate to a spectrum
        weight : np.ndarray, optional
            Dimensions (1, n) for weighting the spectrum of each object

        """

        self.lam = lam
        self.spectra = spectra
        self.x = x
        self.y = y
        self.ref = ref
        if weight is not None:
            self.weight = weight
        else:
            self.weight = np.array([1] * len(x))
        self.lam_res = np.median(lam[1:] - lam[:-1])

        if len(spectra.shape) == 1:
            self.spectra = np.array([spectra])

        self._convert_to_photons()

    def __str__(self):
        return "A photon source object"

    def __getitem__(self, i):
        return (self.x[i], self.y[i],
                self.spectra[self.ref[i], :] * self.weight[i])

    def __mul__(self, x):
        newsrc = deepcopy(self)
        if isinstance(x, (TransmissionCurve, EmissionCurve,
                          UnityCurve, BlackbodyCurve)):
            newsrc._apply_transmission_curve(x)
        else:
            newsrc.spectra *= x
        return newsrc

    def __add__(self, x):
        newsrc = deepcopy(self)
        if isinstance(x, Source):
            if self.units != x.units:
                raise ValueError("units are not compatible: " +
                                 str(self.units) + ", " + str(x.units))

            newsrc.lam = self.lam
            newsrc.spectra = list(self.spectra)
            # Resample new spectra to wavelength grid of self
            for spec in x.spectra:
                tmp = np.interp(self.lam, x.lam, spec)
                scale_factor = np.sum(spec) / np.sum(tmp)
                newsrc.spectra += [tmp * scale_factor]
            newsrc.spectra = np.asarray(newsrc.spectra)
            newsrc.spectra_orig = newsrc.spectra
            newsrc.x = np.array((list(self.x) + list(x.x)))
            newsrc.y = np.array((list(self.y) + list(x.y)))
            newsrc.ref = np.array((list(self.ref) +
                                   list(x.ref + self.spectra.shape[0])))
            newsrc.weight = np.array((list(self.weight) + list(x.weight)))

            newsrc.x_orig = deepcopy(newsrc.x)
            newsrc.y_orig = deepcopy(newsrc.y)

        else:
            newsrc.spectra += x

        newsrc.info["object"] = "combined"

        return newsrc

    def __sub__(self, x):
        newsrc = deepcopy(self)
        newsrc.spectra -= x
        return newsrc

    def __rmul__(self, x):
        return self.__mul__(x)

    def __radd__(self, x):
        return self.__add__(x)

    def __rsub__(self, x):
        return self.__mul__(-1) + x

    def __imul__(self, x):
        return self.__mul__(x)

    def __iadd__(self, x):
        return self.__add__(x)

    def __isub__(self, x):
        return self.__sub__(x)


def load(filename):
    """Load :class:'Source' object from filename"""
    return Source.load(filename)