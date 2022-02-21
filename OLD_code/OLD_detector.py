"""
A description of the Chip noise properties and their positions on the Detector

Module Summary
--------------
This module holds three classes: ``Detector``, ``Chip`` and ``HXRGNoise``.

``Chip``
Everything to do with photons and electrons happens in the ``Chip`` class. Each
``Chip`` is initialised with a position relative to the centre of the detector
array, a size [in pixels] and a resolution [in arcsec]. Photons fall onto the
``Chip`` s and are read out together with the read noise characteristics of the
``Chip``.

``Detector``
The ``Detector`` holds the information on where the ``Chip`` s are placed on the
focal plane. Focal plane coordinates are in [arcsec]. These coordinates are
either read in from a default file or determined by the user. The ``Detector``
object is an intermediary - it only passes information on the photons to the
``Chip`` s. It is mainly a convenience class so that the user can read out all
``Chip`` s at the same time.


Classes
-------
Detector
    builds an array of ``Chip`` py_objects based on a ``UserCommands`` object
Chip
    converts incoming photons into ADUs and adds in read-out noise

See Also
--------
OpticalTrain, Source

Notes
-----

References
----------
[1] Bernhard Rauscher's HxRG Noise Generator script


Examples
--------
The :class:`.Detector` can be used in stand-alone mode. In this case it outputs
only the noise that a sealed-off detector would generate:

    >>> import scopesim
    >>> fpa = scopesim.Detector(scopesim.UserCommands())
    >>> fpa.read_out(output=True, chips=[0])


The :class:`.Detector` is more useful if we combine it with a :class:`.Source`
object and an :class:`.OpticalTrain`. Here we create a :class:`.Source` object
for an open cluster in the LMC and pass the photons arriving from it through the
E-ELT and MICADO. The photons are then cast onto the detector array. Each
:class:`.Chip` converts the photons to ADUs and adds the resulting image to an
Astropy ``HDUList``. The ``HDUList`` is then written to disk.

    >>> # Create a set of tests_commands, optical train and detector
    >>>
    >>> import scopesim
    >>> cmds = scopesim.UserCommands()
    >>> opt_train = scopesim.OpticalTrain(cmds)
    >>> fpa = scopesim.Detector(cmds)
    >>>
    >>> # Pass photons from a 10^4 Msun open cluster in the LMC through to the detector
    >>>
    >>> st = sim.source.source_1E4_Msun_cluster()
    >>> st.apply_optical_train(opt_train, fpa)
    >>>
    >>># Read out the detector array to a FITS file
    >>>
    >>> fpa.read_out(filename="my_raw_image.fits")

"""

############################
#       TODO
# - update open, write - remove references to self.params


import os
import sys
from datetime import datetime

import logging
from copy import deepcopy

import multiprocessing as mp

import numpy as np
import numpy.ma as ma
import scipy.ndimage.interpolation as spi
#from scipy.ndimage.interpolation import zoom

from astropy.io import fits
from astropy.io import ascii as ioascii  # ascii redefines builtin
from astropy.wcs import WCS

#from astropy.stats.funcs import median_absolute_deviation as mad

from scopesim.utils import find_file, airmass2zendist

from OLD_code import OLD_simcado_user_commands, OLD_spectral as sc
from scopesim.detector.nghxrg import HXRGNoise

__all__ = ["Detector", "Chip", "open", "plot_detector", "plot_detector_layout",
           "make_noise_cube", "install_noise_cube"]


################################################################################
#                              Detector Objects                                #
################################################################################

class Detector(object):
    """
    Generate a series of :class:`.Chip` py_objects for a focal plane array


    Summary
    -------
    The :class:`.Detector` is a holder for the series of :class:`.Chip`
    py_objects which make up the detector array. The main advantage of the
    :class:`.Detector` object is that the  user can read out all chips in the
    whole detector array at once. A :class:`.Detector` is a parameter in the
    :meth:`.Source.apply_optical_train()` method.


    Parameters
    ----------
    cmds : UserCommands
        Commands for how to model the Detector
    small_fov : bool, optional
        Default is True. Uses a single 1024x1024 window at the centre of the FoV

    Attributes
    ----------
    cmds : UserCommands
        tests_commands for modelling the detector layout and exposures
    layout : astropy.table.Table
        table of positions and sizes of the chips on the focal plane
    chips : list
        a list of the ``Chips`` which make up the detector array
    oversample : int
        factor between the internal angular resolution and the pixel FOV
    fpa_res : float
        [mas] field of view of a single pixel
    exptime : float
        [s] exposure time of a single DIT
    tro : float
        [s] time between consecutive non-destructive readouts in up-the-ramp mode
    ndit : int
        number of exposures (DITs)


    Methods
    -------
    read_out()
        for reading out the detector array into a FITS file
    open()
        not yet implemented
    write()
        not yet implemented
        Save the Detector object into a FITS file

    .. todo::

        Open should be moved into a general function for OLD_detector.py which
        returns a :class:`.Detector` object after reading in a saved detector
        file


    See Also
    --------
    :class:`.Chip`, :class:`.Source`,
    :class:`.OpticalTrain`, :class:`.UserCommands`

    Examples
    --------
    Create a :class:`Detector` object
    ::

        >>> import scopesim as sim
        >>> my_cmds = sim.UserCommands()
        >>> my_detector = sim.Detector(my_cmds)


    Read out only the first :class:`.Chip`
    ::

        >>> my_detector.readout(filename=image.fits, chips=[0])


    """


    def __init__(self, cmds, small_fov=True):
        # 1. Read in the chip layout
        # 2. Read in the flat field file
        # 3. Generate chip py_objects
        # 4. Check if a noise file has been given
            # if not, generate new noise files
            # else: read in the noise file
            # if the noise file has many extensions, choose several random
            #    extensions

        self.cmds = cmds

        if small_fov:
            print("Safety switch is on - Detector(..., small_fov='True')")
            self.layout = ioascii.read(
                """#  id    x_cen    y_cen   x_len   y_len   pixsize  angle    gain
                   #           mm       mm   pixel   pixel        mm    deg  e-/ADU
                       0        0        0    1024    1024     0.015      0.      1.""")
        else:
            try:
                self.layout = ioascii.read(self.cmds["FPA_CHIP_LAYOUT"])
            except:
                raise FileNotFoundError(self.cmds["FPA_CHIP_LAYOUT"] +
                                        " (FPA_CHIP_LAYOUT) cannot be read")


        if self.cmds["INST_FLAT_FIELD"] is not None and \
           os.path.exists(self.cmds["INST_FLAT_FIELD"]):
            hdu_flat_field = fits.open(self.cmds["INST_FLAT_FIELD"])
            n_exts = len(hdu_flat_field.info(output=False))
            if n_exts == len(self.layout["id"]):
                j = 0
            elif n_exts > len(self.layout["id"]):
                j = 1
            elif n_exts == 1:
                hdu_flat_field = [hdu_flat_field[0]] * len(self.layout["id"])
            else:
                raise ValueError("What's going on with INST_FLAT_FIELD?")
        else:
            j = 0
            hdu_flat_field = [None] * len(self.layout["id"])


        self.chips = [Chip(self.layout["x_cen"][i], self.layout["y_cen"][i],
                           self.layout["x_len"][i], self.layout["y_len"][i],
                           self.cmds["SIM_PIXEL_SCALE"],
                           self.layout["pixsize"][i],
                           self.layout["angle"][i],
                           self.layout["gain"][i],
                           (self.cmds["OBS_RA"], self.cmds["OBS_DEC"]),
                           0.,
                           self.layout["id"][i],
                           hdu_flat_field[i+j])
                      for i in range(len(self.layout["x_cen"]))]

        self.oversample = self.cmds["SIM_OVERSAMPLING"]
        self.fpa_res = self.cmds["SIM_PIXEL_SCALE"]
        self.exptime = self.cmds["OBS_EXPTIME"]
        self.ndit    = self.cmds["OBS_NDIT"]
        self._n_ph_atmo   = 0
        self._n_ph_mirror = 0
        self._n_ph_ao     = 0
        self.array = None        # defined in method


    def read_out(self, filename=None, to_disk=False, chips=None,
                 read_out_type="superfast", **kwargs):
        """
        Simulate the read-out process of the detector array

        Based on the parameters set in the ``UserCommands`` object, the detector
        will read out the images stored on the ``Chips`` according to the
        specified read-out scheme, i.e. Fowler, up-the-ramp, single read, etc.

        Parameters
        ----------
        filename : str
            where the file is to be saved. If ``None`` and ``to_disk`` is true,
            the output file is called "output.fits". Default is ``None``

        to_disk : bool
            a flag for where the output should go.
            If ``filename`` is given or if ``to_disk=True``,
            the ``Chip`` images will be written to a  `.fits`` file on disk.
            If no `filename`` is specified, the output will be called "output.fits".

        chips : int, array-like, optional
            The chip or chips to be read out, based on the detector_layout.dat
            file. Default is the first ``Chip`` specified in the list, i.e. [0].

        read_out_type : str, optional
            The name of the algorithm used to read out the chips:
            - "superfast"
            - "non_destructive"
            - "up_the_ramp"

        Returns
        -------
        astropy.io.fits.HDUList

        Keyword Arguments (**kwargs)
        ----------------------------
        **kwargs are used to update the ``UserCommands`` object that controls
        the ``Detector``. Therefore any dictionary keywords can be passed in the
        form of a dictionary, i.e. {"OBS_EXPTIME" : 60, "OBS_OUTPUT_DIR" : "./"}

        """

        #removed kwargs
        self.cmds.update(kwargs)

        if filename is not None:
            to_disk = True

        if filename is None and to_disk is True:
            if self.cmds["OBS_OUTPUT_DIR"] is None:
                self.cmds["OBS_OUTPUT_DIR"] = "./output.fits"
            filename = self.cmds["OBS_OUTPUT_DIR"]

        if chips is not None:
            if np.isscalar(chips):
                ro_chips = [chips]
            else:
                ro_chips = chips
        elif chips is None:
            ro_chips = np.arange(len(self.chips))
        else:
            raise ValueError("Something wrong with ``chips``")

        # Time stamp for FITS header
        #creation_date = datetime.now().isoformat(timespec='seconds')
        # timespec="seconds" throws an error on some python versions
        creation_date = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

        hdulist = fits.HDUList()

        # Create primary header unit for multi-extension files
        if len(ro_chips) > 1:
            primary_hdu = fits.PrimaryHDU()

            primary_hdu.header['DATE'] = creation_date


            for key in self.cmds.cmds:
                val = self.cmds.cmds[key]

                if isinstance(val, (sc.TransmissionCurve, sc.EmissionCurve,
                                    sc.UnityCurve, sc.BlackbodyCurve)):
                    val = val.params["filename"]

                if isinstance(val, str) and len(val) > 35:
                    val = "... " + val[-35:]

                try:
                    primary_hdu.header["HIERARCH "+key] = val
                except NameError:   # any other exceptions possible?
                    pass
            hdulist.append(primary_hdu)

        # Save the detector image(s)
        for i in ro_chips:
            ######
            # Put in a catch here so that only the chips specified in "chips"
            # are read out
            ######
            print("Reading out chip", self.chips[i].id, "using",
                  read_out_type)

            array = self.chips[i].read_out(self.cmds,
                                           read_out_type=read_out_type)

            ## TODO: transpose is just a hack - need to make sure
            ##       x and y are handled correctly throughout ScopeSim
            thishdu = fits.ImageHDU(array.T)

            thishdu.header["EXTNAME"] = ("CHIP_{:02d}".format(self.chips[i].id),
                                         "Chip ID")

            thishdu.header["CHIP_ID"] = (self.chips[i].id, "Chip ID")
            thishdu.header['DATE'] = creation_date

            # Primary WCS for sky coordinates
            thishdu.header.extend(self.chips[i].wcs.to_header())

            # Secondary WCS for focal plane coordinates
            try:
                thishdu.header.extend(self.chips[i].wcs_fp.to_header(key='A'))
            except AttributeError:
                print("No WCS_FP!")
                pass

            thishdu.header["BUNIT"] = ("ADU", "")
            thishdu.header["EXPTIME"] = (self.exptime, "[s] Exposure time")
            thishdu.header["NDIT"] = (self.ndit, "Number of exposures")
            #thishdu.header["TRO"] = (self.tro,
            #                         "[s] Time between non-destructive readouts")
            thishdu.header["GAIN"] = (self.chips[i].gain, "[e-/ADU]")
            thishdu.header["AIRMASS"] = (self.cmds["ATMO_AIRMASS"], "")
            thishdu.header["ZD"] = \
                (airmass2zendist(self.cmds["ATMO_AIRMASS"]), "[deg]")


            for key in self.cmds.cmds:
                val = self.cmds.cmds[key]
                if isinstance(val, str):
                    if len(val) > 35:
                        val = "... " + val[-35:]
                try:
                    thishdu.header["HIERARCH "+key] = val
                except NameError:   # any other exceptions possible?
                    pass
                except ValueError:
                    logging.warning("ValueError - Couldn't add keyword: "+key)
            hdulist.append(thishdu)

        if to_disk:
            hdulist.writeto(filename, clobber=True, checksum=True)

        return hdulist


    def write(self, filename=None, **kwargs):
        """
        Write a ``Detector`` object out to a FITS file


        Summary
        -------
        Writes the important information contained in a ``Detector`` object into
        FITS file for later use. The main information written out includes: the
        layout of the detector chips, any pixel maps associated with the
        detector chips, a linearity curve and a QE curve for the chips.


        Parameters
        ----------
        filename : str, optional
            path to the FITS file where the ``Detector`` object is stored. If
            ``filename=None`` (by default), the file written is ``./detector.fits``



        Keyword Arguments (**kwargs)
        ----------------------------

        Examples
        --------
        """

        self.params.update(kwargs)
        if filename is None:
            filename = self.params["OBS_OUTPUT_DIR"]
        if filename is None:
            raise ValueError("No output path was specified. " + \
                             "Use either filename= or HXRG_OUTPUT_PATH=")

        hdu = fits.PrimaryHDU(self.array)
        hdu.header["PIX_RES"] = self.pix_res
        hdu.header["EXPTIME"] = self.exptime
        hdu.header["GAIN"]    = self.params["FPA_GAIN"]
        hdu.header["SIMCADO"] = "FPA_NOISE"

        try:
            hdu.writeto(filename, clobber=True, checksum=True)
        except OSError:
            logging.warning(filename+" exists and is busy. OS won't let me write")


################################################################################
#                              Chip Objects                                    #
################################################################################


class Chip(object):
    """
    Holds the "image" as seen by a single chip in the focal plane


    Summary
    -------
    The ``Chip`` object contains information on where it is located in the focal
    plane array. The method ``<Source>.apply_optical_train()`` passes an image of
    the on-sky object to each ``Chip``. This image is resampled to the ``Chip``
    pixel scale. Each ``Chip`` holds the "ideal" image as an array of expectation
    values for the level of photons arriving during an EXPTIME. The ``Chip`` then
    adds detector noise and other characteristics to the image when
    <Detector>.readout() is called.


    Parameters
    ----------
    x_cen, y_cen : float
        [micron] the coordinates of the centre of the chip relative to the
        centre of the focal plane
    x_len, y_len : int
        the number of pixels per dimension
    pix_res : float
        [arcsec] the field of view per pixel
    id : int
        an identification number for the chip (assuming they are not correctly
        ordered)
    flat_field : np.ndarray
        a 2D array holding the flat fielding effects for the chip

    Attributes
    ----------
    x_cen, y_cen : float
        [arcsec] the coordinates of the centre of the chip relative to the
        centre of the focal plane
    naxis1, naxis2 : int
        the number of pixels per dimension
    pix_res : float
        [arcsec] the field of view per pixel
    chipid : int, optional
        the id of the chip relative to the others on the detector array. Default is ``None``
    dx, dy : float
        [arcsec] half of the field of view of each chip
    x_min, x_max, y_min, y_max : float
        [arcsec] the borders of the chip relative to the centre of the focal plane
    array : np.ndarray
        an array for holding the signal registered by the ``Chip``


    Methods
    -------
    add_signal(signal)
        adds signal to ``.array``. The signal should be the same dimensions as
        ``Chip.array``
    add_uniform_background(emission, lam_min, lam_max, output=False)
        adds a constant to the signal in ``.array``. The background level is found
        by integrating the ``emission`` curve between ``lam_min`` and ``lam_max``.
        If output is set to ``True``, an image with the same dimensions as
        ``.array`` scaled to the background flux is returned.
    apply_pixel_map(pixel_map_path=None, dead_pix=None, max_well_depth=1E5)
        applies a mask to ``.array`` representing the position of the current
        "hot" and "dead" pixels / lines
    reset()
        resets the signal on the ``Chip`` to zero. In future releases, an
        implementation of the persistence characteristics of the detector will
        go here.


    Raises
    ------

    See Also
    --------
    Detector, Source, UserCommands, OpticalTrain

    Examples
    --------
    """

    def __init__(self, x_cen, y_cen, x_len, y_len, pix_res, pixsize=15,
                 angle=0, gain=1, obs_coords=[0, 0], fieldangle=0,
                 chipid=None, flat_field=None):

        # permanent rotation of chip
        cangle = np.cos(np.deg2rad(angle))
        sangle = np.sin(np.deg2rad(angle))

        # offset of chip centre from field centre in um
        xoff = x_cen * cangle - y_cen * sangle
        yoff = x_cen * sangle + y_cen * cangle

        # Primary WCS to transform pixel coordinates to sky coordinates
        thewcs = WCS(naxis=2)
        thewcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        thewcs.wcs.name = "PIX2SKY"
        thewcs.wcs.crval = obs_coords
        thewcs.wcs.crpix = [(x_len + 1) / 2 - xoff / pixsize,
                            (y_len + 1) / 2 - yoff / pixsize]
        thewcs.wcs.cdelt = [pix_res / 3600, pix_res / 3600]  # arcsec to deg
        thewcs.wcs.cunit = ['deg', 'deg']
        thewcs.wcs.pc = np.array([[cangle, sangle], [-sangle, cangle]])

        fieldangle = np.deg2rad(fieldangle)
        fieldrot = np.matrix([[np.cos(fieldangle), np.sin(fieldangle)],
                              [-np.sin(fieldangle), np.cos(fieldangle)]])
        thewcs.wcs.pc = fieldrot * np.asmatrix(thewcs.wcs.pc)

        self.wcs = thewcs

        # Second WCS to transform pixel coordinates to focal-plane coordinates
        wcs_fp = WCS(naxis=2)
        wcs_fp.wcs.ctype = ['LINEAR', 'LINEAR']
        wcs_fp.wcs.name = "PIX2FP"
        wcs_fp.wcs.crval = [0., 0.]
        wcs_fp.wcs.crpix = [(x_len + 1) / 2 - xoff / pixsize,
                            (y_len + 1) / 2 - yoff / pixsize]
        wcs_fp.wcs.cdelt = [pixsize, pixsize]
        wcs_fp.wcs.cunit = ['mm', 'mm']
        wcs_fp.wcs.pc = np.array([[cangle, sangle], [-sangle, cangle]])

        self.wcs_fp = wcs_fp

        # attributes in um
        self.xcen_um = x_cen
        self.ycen_um = y_cen
        dx_um = (x_len // 2) * pixsize
        dy_um = (y_len // 2) * pixsize
        self.xmin_um = self.xcen_um - dx_um
        self.xmax_um = self.xcen_um + dx_um
        self.ymin_um = self.ycen_um - dy_um
        self.ymax_um = self.ycen_um + dy_um

        # for backwards compatibility: attributes in arcsec
        self.x_cen  = x_cen / pixsize * pix_res
        self.y_cen  = y_cen / pixsize * pix_res
        self.naxis1 = x_len
        self.naxis2 = y_len
        self.pix_res = pix_res
        self.gain   = gain
        self.id     = chipid

        dx = (x_len // 2) * pix_res
        dy = (y_len // 2) * pix_res
        self.x_min = self.x_cen - dx
        self.x_max = self.x_cen + dx
        self.y_min = self.y_cen - dy
        self.y_max = self.y_cen + dy

        self.array = None

        self.ndit    = 0
        self.dit     = 0
        self.dark    = 0
        self.min_dit = 0

        if flat_field is not None:
            if isinstance(flat_field, (fits.ImageHDU, fits.PrimaryHDU)):
                flat_field = flat_field.data
            wf, hf = flat_field.shape
            zx, zy = self.naxis1 / wf, self.naxis2 / hf
            if wf == 1. and hf == 1.:
                self.flat_field = flat_field
            else:
                self.flat_field = spi.zoom(flat_field, (zx, zy), order=1)
            self.flat_field = self.flat_field.T
        else:
            self.flat_field = None

    def add_signal(self, signal):
        """
        Add a 2D array of photon signal to the Chip

        Summary
        -------
        Add some signal photons to the detector array. Input units are expected
        to be [ph/s/pixel]

        Parameters
        ----------
        signal : np.ndarray
            [ph/pixel/s] photon signal. ``signal`` should have the same dimensions
            as the ``array``

        Returns
        -------
        None

        """
        if isinstance(signal, np.ndarray):
            if signal.shape[0] == self.naxis1 and \
               signal.shape[1] == self.naxis2:
                if self.array is None:
                    self.array = signal
                else:
                    self.array += signal
            else:
                raise ValueError(str(signal.shape) + " != " + \
                                 str((self.naxis1, self.naxis2)))
        elif not hasattr(signal, "__len__"):
            self.array += signal

    def add_uniform_background(self, emission, lam_min, lam_max, output=False):
        """
        Add a uniform background

        Summary
        -------
        Take an EmissionCurve and some wavelength boundaries, lam_min and lam_max,
        and sum up the photons in between. Add those to the source array.

        Parameters
        ----------
        - emission_curve: EmissionCurve object with background emission photons
        - lam_min, lam_max: the wavelength limits

        Optional keywords:
        - output: [False, True] if output is True, the BG emission array is
                  returned

        Output is in [ph/s/pixel].
        """

        if isinstance(emission, sc.EmissionCurve):
            bg_photons = emission.photons_in_range(lam_min, lam_max)
        elif isinstance(emission, (float, int)):
            bg_photons = emission
        else:
            bg_photons = 0
            logging.warning("type(emission) invalid. No background added")

        if output is True:
            return bg_photons * np.ones(self.array.shape, dtype=np.float32)
        else:
            self.array += bg_photons


    def apply_pixel_map(self, pixel_map_path=None, dead_pix=None,
                        max_well_depth=1E5):
        """
        Adds "hot" and "dead" pixels to the array

        Summary
        -------
        applies a mask to ``.array`` representing the positions of the current
        "hot" and "dead" pixels / lines. The method either reads in a FITS file
        with locations of these pixels, or generates a series of random
        coordinates and random weights for the pixels.

        Parameters
        ----------
        pixel_map_path : str
            path to the FITS file. Default is None
        dead_pix : int
            [%] the percentage of dead or hot pixels on the chip - only used if
            ``pixel_map_path = None``. Default is ``None``.
        max_well_depth : 1E5


        Returns
        -------
        None

        """

        try:
            pixel_map = fits.getdata(pixel_map_path)
            if self.array.shape != pixel_map.shape:
                raise ValueError("pixel_map.shape != detector_array.shape")
            self.array += pixel_map * max_well_depth
        except ValueError:
            if dead_pix is not None:
                n = int(self.naxis1 * self.naxis2 * dead_pix / 100)
                x = np.random.randint(self.naxis1, size=n)
                y = np.random.randint(self.naxis2, size=n)
                z = np.random.random(n)
                self.array[x, y] += z * max_well_depth
            else:
                raise ValueError("Couldn't apply pixel_map")


    def reset(self):
        self.array = None


    def read_out(self, cmds, read_out_type="superfast"):
        """
        Read out the detector array

        Parameters
        ----------
        cmds : scopesim.UserCommands
            Commands for how to read out the chip

        Returns
        -------
        out_array : np.ndarray
            image of the chip read out

        """

        # set up the read out
        self.dit      = cmds["OBS_EXPTIME"] / cmds["OBS_NDIT"]
        self.ndit     = int(cmds["OBS_NDIT"])
        self.dark     = cmds["FPA_DARK_MEDIAN"]
        self.min_dit  = cmds["FPA_PIXEL_READ_TIME"] * \
                        (self.naxis1 * self.naxis1 / cmds["HXRG_NUM_OUTPUTS"])

        if self.array is None:
            self.array = np.zeros((self.naxis2, self.naxis1), dtype=np.float32)

        # At this point, the only negatives come from the convolution.
        # Remove them for the Poisson process
        self.array[self.array < 0] = 0
        out_array = np.zeros(self.array.shape, dtype=np.float32)

        ######## Multiply by Exptime
        # the different read out modes
        if read_out_type.lower() == "superfast":
            out_array = self._read_out_superfast(cmds, self.dit, self.ndit)
        elif read_out_type.lower() == "non_destructive":
            out_array = self._read_out_non_destructive(cmds, self.dit, self.ndit)
        elif read_out_type.lower() == "up_the_ramp":
            out_array = self._read_out_up_the_ramp(cmds, self.dit, max_byte=2**30)
        else:
            raise ValueError("``read_out_type`` not readable")

        ######## Flat fielding
        if self.flat_field is not None:
            out_array *= self.flat_field

        ######## Remove a constant BG level
        if cmds["OBS_REMOVE_CONST_BG"].lower() == "yes":
            bg_val = np.median(out_array)
            out_array -= bg_val

        return out_array

    def _read_out_non_destructive(self, cmds, dit, ndit):
        """
        Read out NDIT times non-destructively according to FPA_READ_OUT_SCHEME

        Parameters
        ----------
        cmds : UserCommands

        Returns
        -------
        out_array : np.ndarray

        """
        lin_curve = cmds["FPA_LINEARITY_CURVE"]
        ro_times  = self._get_readout_times(scheme=cmds["FPA_READ_OUT_SCHEME"])
        out_array = np.zeros(self.array.shape, dtype=np.float32)

        for n in range(self.ndit):
            ro_cube = []
            for t in ro_times:

                signal = self._read_out_poisson((self.array + self.dark),
                                                dit=t, ndit=1)
                if lin_curve is not None:
                    signal, lin_curve = self._apply_linearity(signal, lin_curve,
                                                              return_curve=True)
                read_noise = self._read_noise_frame(cmds)[0]
                ro_cube += [signal + read_noise]

            ro_cube = np.array(ro_cube)
            mask = ro_cube > cmds["FPA_FULL_WELL_DEPTH"]
            masked_ro_cube = ma.array(ro_cube, mask=mask)
            masked_ro_cube[1:,:,:] -= masked_ro_cube[0,:,:]
            ro_times.resize((len(ro_times),1,1))
            masked_ro_cube[1:,:,:] /= ro_times[1:,:,:]
            av_ro_cube = np.average(masked_ro_cube[1:,:,:], axis=0)

            out_array += av_ro_cube * dit

        out_array /= self.gain

        return out_array.astype(np.float32)

    # TODO: What to do if dit = min_dit (single read)?
    # TODO: Make breaking up into memory chunks more flexible?
    def _read_out_up_the_ramp(self, cmds, dit, max_byte=2**30):
        """
        Test readout onto a detector using cube model

        Parameters
        ----------
        image :
            a 2D image to be mapped onto the detector. Units are [ph/s/pixel]
        dit :
            integration time [s]

        Optional Parameters
        -------------------
        max_byte :
            the largest possible chunk of memory that can be used for computing
            the sampling slope

        This function builds an intermediate cube of dimensions (nx, ny, nro) with a
        layer for  each non-destructive read.

        Output is given in [ph/pixel].
        """

        print("NOTE - 'up the ramp' readout only reads a single DIT")

        image = self.array
        tro   = self.min_dit

        nx, ny = image.shape

        nro = np.int(dit / tro)
        tpts = (1 + np.arange(nro)) * tro

        img_byte = image.nbytes
        pix_byte = img_byte / (nx * ny)

        max_pix = max_byte / pix_byte

        #cube_megabyte = img_byte * nro / 2**20
        #print("Full cube  has {0:.1f} Megabytes".format(cube_megabyte))

        ny_cut = np.int(max_pix / (nx * nro))
        if ny_cut >= ny:
            ny_cut = ny
        #print("Cut image to ny={0:d} rows".format(ny_cut))

        slope = np.zeros(image.shape)
        cube = np.zeros((nx, ny_cut, nro))

        y1 = 0
        while y1 < ny:

            y2 = y1 + ny_cut
            if y2 > ny:
                y2 = ny
                ny_cut = ny - y1
                try:
                    del cube
                    cube = np.zeros((nx, ny_cut, nro))
                except NameError:
                    pass

            ## Fill the cube with Poisson realization, individual reads
            for i in range(nro):
                cube[:, :, i] = np.random.poisson(image[:, y1:y2] * tro)

            ## Build the ramp
            sumcube = cube.cumsum(axis=2)

            ## determine the slope using explicit formula calculated over cube
            Sx = tpts.sum()
            Sxx = (tpts * tpts).sum()
            Sy = np.sum(sumcube, axis=2)
            Sxy = np.sum(sumcube * tpts, axis=2)

            slope[:, y1:y2] = (nro * Sxy - Sx * Sy) / (nro * Sxx - Sx * Sx)

            ## Move to next slice
            y1 = y2

        # return values are [ph/pixel]
        return slope * dit


    def _read_out_superfast(self, cmds, dit, ndit):
        """
        Superfast read-out
        """

        signal = self._read_out_poisson(self.array, dit, ndit)

        # apply the linearity curve
        lin_curve = cmds["FPA_LINEARITY_CURVE"]
        if lin_curve is not None:
            signal, lin_curve = self._apply_linearity(signal,
                                                      lin_curve,
                                                      return_curve=True)

        # superfast hack to get an approximation of the readout noise
        # in the image
        ro = self._read_noise_frame(cmds, n_frames=1) * np.sqrt(ndit)

        ########## Could work, but it's too slow for ndit > 10 ##############
        # add 1 to the ndits, because there will always be a readout at
        # the start
        #ro_frames = self._read_noise_frame(cmds, n_frames=max(ndit,2))
        #ro = np.sum(ro_frames, axis=0)

        out_array = signal + ro
        out_array /= self.gain

        return out_array


    def _read_out_poisson(self, image, dit, ndit):
        """
        Apply a poisson distribution to the image

        Parameters
        ----------
        image : np.ndarray
            image to be poissonified
        dit : float
            [s] length of exposure
        ndit : int
            [#] number of exposures

        Returns
        -------
        im_st : np.ndarray
            average of ndit exposures of length dit

        """
        image2 = image * dit
        ## does not seem to be necessary in numpy version 1.12.1 any more
#        image2[image2 > 2.14E9] = 2.14E9

        im_st = np.zeros(np.shape(image))
        for _ in range(ndit):
            im_st += np.random.poisson(image2)

        return im_st.astype(np.float32)


    def _read_noise_frame(self, cmds, n_frames=1):
        """
        Read in read-out-noise from the FITS file specified by FPA_NOISE_PATH

        If cmds["FPA_NOISE_PATH"] == "gen", all info about the size and number
        of layers must be in the :class:`.UserCommands` object.
        I.e. HXRG_NUM_NDRO

        Parameters
        ----------
        cmds : UserCommands
        n_frames : int
            The number of frames needed

        Returns
        -------
        noise_cube : np.ndarray
            if n_frames == 1: shape = (naxis1, naxis2)
            else: shape = = (naxis1, naxis2, n_frames)

        """

        if cmds["FPA_USE_NOISE"].lower() == "no":
            return np.zeros((self.naxis1, self.naxis2))

        if "gen" in cmds["FPA_NOISE_PATH"].lower():
            if cmds["HXRG_OUTPUT_PATH"] is not None:
                generate_hxrg_noise(cmds)
                tmp = fits.getdata(cmds["HXRG_OUTPUT_PATH"])
                return tmp[:self.naxis1, :self.naxis2]
            else:
                noise_cube = generate_hxrg_noise(cmds)

        elif cmds["FPA_NOISE_PATH"] is not None:
            n = len(fits.info(cmds["FPA_NOISE_PATH"], False))
            layer = np.random.randint(low=0, high=n, size=n_frames)
            tmp = [fits.getdata(cmds["FPA_NOISE_PATH"], i) for i in layer]
            noise_cube = np.array([im[:self.naxis1, :self.naxis2] for im in tmp])

        else:
            noise_cube = np.zeros((self.naxis1, self.naxis2, n_frames))

        if n_frames == 1:
            return noise_cube[0,:,:]
        else:
            return noise_cube


    def _get_readout_times(self, scheme="double_corr"):
        """
        Expect that scheme = cmds["FPA_READ_OUT_SCHEME"]

        """

        if "double_corr" in scheme:
            scheme = [0, self.dit]
        elif "up" in scheme:
            scheme = np.arange(0, self.dit, self.min_dit + 1E-3)
        elif "fowl" in scheme:
            fowl_pair = np.min((4, int(self.dit / self.min_dit) // 2))
            scheme = np.arange(0, self.dit, self.min_dit + 1E-3).tolist()
            scheme = scheme[:fowl_pair] + scheme[-fowl_pair:]


        if isinstance(scheme, str):
            if os.path.exists(scheme):
                data = ioascii.read(scheme)
                times = data[data.colnames[0]]
        elif isinstance(scheme, (list, tuple, np.ndarray)):
            times = np.array(scheme)

        return times.astype(np.float32)


    def _apply_linearity(self, in_array, curve, return_curve=False):
        """
        Apply a linearity curve to a 2D array

        Parameters
        ----------
        in_array : array
            the array to be linearized
        curve : string, astropy.Table
            if string, it is assumed to be a filename
            if Table,  it is assumed to be from a previous run
        return_curve : bool
            whether to return the linearity curve

        Returns
        -------
        out_array, data : array, Table
            if return_curve == True
        out_array : array
            if return_curve == False

        """

        from astropy.table import Table

        if isinstance(curve, str):
            if os.path.exists(curve):
                data = ioascii.read(curve)
            else:
                raise ValueError("file doesn't exist: "+curve)
        elif isinstance(curve, Table):
            data = curve

        real_cts = data[data.colnames[0]].data.astype(np.float32)
        measured_cts = data[data.colnames[1]].data.astype(np.float32)

        out_array = np.interp(in_array.flatten(),
                              real_cts, measured_cts).reshape(in_array.shape)
        out_array = out_array.astype(np.float32)

        if return_curve:
            return out_array, data
        else:
            return out_array




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



# TODO this ought to be renamed (redefined-builtin)
def open(self, filename):
    """
    Opens a saved ``Detector`` file.


    Summary
    -------
    ** Not yet implemented **
    ** Should be moved outside of ``Detector`` and called with
    ``detector.open()`` **

    Detector py_objects can be saved to FITS file and read back in for later
    simulations.

    Parameters
    ----------
    filename : str
        path to the FITS file where the ``Detector`` object is stored

    Returns
    -------
    ``scopesim.Detector`` object

    Examples
    --------
    """

    if not os.path.exists(filename):
        raise FileNotFoundError(filename + " doesn't exist")

    with fits.open(filename) as fp1:
        self.params.update(fp1[0].header)
        self.array = fp1[0].data

    raise ValueError("Function not finished")


#def plot_detector_layout(detector, clr="g"):
#    """Plot the detector layout. NOT FINISHED """
#    try:
#        import matplotlib.pyplot as plt
#    except ImportError:
#        raise ImportError("matplotlib can't be found")
#
#    for i, chip in enumerate(detector.chips):
#        plt.plot((chip.x_min, chip.x_max), (chip.y_min, chip.y_min), clr)
#        plt.plot((chip.x_min, chip.x_max), (chip.y_max, chip.y_max), clr)
#        plt.plot((chip.x_min, chip.x_min), (chip.y_min, chip.y_max), clr)
#        plt.plot((chip.x_max, chip.x_max), (chip.y_min, chip.y_max), clr)
#        plt.text(chip.x_cen, chip.y_cen, chip.id, fontsize=14)
#        plt.xlabel("Distance [arcsec]", fontsize=14)
#        plt.ylabel("Distance [arcsec]", fontsize=14)

def plot_detector_layout(detector, plane="sky", clr='g-', plot_origin=False):
    """Plot the detector layout"""

    from matplotlib import pyplot as plt
    npts = 101
    for chip in detector.chips:

        if plane == 'sky':
            thewcs = chip.wcs
            scale = 3600.
            xlabel = 'RA offset (arcsec)'
            ylabel = 'DE offset (arcsec)'
        elif plane == 'fpa':
            thewcs = chip.wcs_fp
            scale = 1.
            xlabel = 'x (mm)'
            ylabel = 'y (mm)'

        xrange = np.linspace(1, chip.naxis1, npts)
        yrange = np.linspace(1, chip.naxis2, npts)
        xpix = np.concatenate((xrange,
                               np.zeros(npts) + chip.naxis1,
                               xrange[::-1],
                               np.zeros(npts) + 1))
        ypix = np.concatenate((np.zeros(npts) + 1,
                               yrange,
                               np.zeros(npts) + chip.naxis2,
                               yrange[::-1]))

        xworld, yworld = thewcs.all_pix2world(xpix, ypix, 1)
        xworld -= thewcs.wcs.crval[0]
        yworld -= thewcs.wcs.crval[1]
        plt.plot(xworld * scale, yworld * scale, clr)

        if plot_origin:
            x0, y0 = thewcs.all_pix2world(1, 1, 1)
            x0 -= thewcs.wcs.crval[0]
            y0 -= thewcs.wcs.crval[1]
            plt.plot(x0 * scale, y0 * scale, 'r.')

        xcen, ycen = thewcs.all_pix2world(chip.naxis1 / 2, chip.naxis2 / 2, 1)
        xcen -= thewcs.wcs.crval[0]
        ycen -= thewcs.wcs.crval[1]
        plt.text(xcen * scale, ycen * scale, chip.id)

    plt.axes().set_aspect('equal')
    if plane == 'sky':
        plt.gca().invert_xaxis()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


def plot_detector(detector):
    """
    Plot the contents of a detector array

    Parameters
    ----------
    detector : scopesim.Detector
        The detector object to be shown
    """

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.figure(figsize=(15*2, 13.75*2))

    x0 = np.min([i.x_min for i in detector.chips])
    x1 = np.max([i.x_max for i in detector.chips])
    y0 = np.min([i.y_min for i in detector.chips])
    y1 = np.max([i.y_max for i in detector.chips])
    w = (detector.chips[0].x_max - detector.chips[0].x_min)/(x1 - x0)
    h = (detector.chips[0].y_max - detector.chips[0].y_min)/(y1 - y0)

    for chip in detector.chips:
        ax1 = plt.axes([(chip.x_min - x0)/(x1 - x0), (chip.y_min - y0)/(y1 - y0),
                      w, h])
        ax1.set_xticks([])
        ax1.set_yticks([])
        if chip.array is not None:
            ax1.imshow(np.rot90(chip.array - np.min(chip.array)), norm=LogNorm(),
                       cmap="Greys", vmin=1)

    plt.show()




def generate_hxrg_noise(cmds):
    """
    Generate a read noise frame using a UserCommands object

    Create a detector noise array using Bernhard Rauscher's NGHxRG tool

    Parameters
    ----------
    cmds : scopesim.UserCommands

    """

    #if len(kwargs) > 0 and self.verbose: print("updating ",kwargs)
    #self.params.update(kwargs)
    print("Generating a new chip noise array")
    print(mp.current_process())
    # HXRG needs a pca file to run. Work out what a PCA file means!!
    ng_h4rg = HXRGNoise(naxis1=cmds["HXRG_NAXIS1"],
                        naxis2=cmds["HXRG_NAXIS2"],
                        naxis3=cmds["HXRG_NUM_NDRO"],
                        n_out=cmds["HXRG_NUM_OUTPUTS"],
                        nroh=cmds["HXRG_NUM_ROW_OH"],
                        pca0_file=cmds["HXRG_PCA0_FILENAME"],
                        verbose=cmds["SIM_VERBOSE"])

    # Make a noise file
    noise = ng_h4rg.mknoise(o_file=cmds["HXRG_OUTPUT_PATH"],
                            rd_noise=cmds["FPA_READOUT_MEDIAN"],
                            pedestal=cmds["HXRG_PEDESTAL"],
                            c_pink=cmds["HXRG_CORR_PINK"],
                            u_pink=cmds["HXRG_UNCORR_PINK"],
                            acn=cmds["HXRG_ALT_COL_NOISE"])

    return noise


def make_noise_cube(num_layers=25, filename="FPA_noise.fits", multicore=True):
    """
    Create a large noise cube with many separate readout frames.

    Note:
    Each frame takes about 15 seconds to be generated. The default value of
    25 frames will take around six minutes depending on your computer's
    architecture.

    Parameters
    ----------
    num_layers : int, optional
        the number of separate readout frames to be generated. Default is 25.
    filename : str, optional
        The filename for the FITS cube. Default is "FPA_noise.fits"
    multicore : bool, optional
        If you're not using windows, this allows the process to use all
        available cores on your machine to speed up the process. Default is True

    Notes
    -----
    multicore doesn't work - fix it

    """

    if sys.version_info.major < 3:
        print("Sorry, but this only works in Python 3 and above. \
           See the ScopeSim FAQs for work-around options")
        return None


    cmds = OLD_simcado_user_commands.UserCommands()
    cmds["FPA_NOISE_PATH"] = "generate"
    cmds["FPA_CHIP_LAYOUT"] = "default"


    if "Windows" in os.environ.get('OS', ''):
        multicore = False

    if __name__ == "__main__" and multicore:
        pool = mp.Pool(processes=mp.cpu_count()-1)
        frames = pool.map(generate_hxrg_noise, (cmds)*num_layers)
    else:
        frames = [generate_hxrg_noise(cmds) \
                  for i in range(num_layers)]

    hdu = fits.HDUList([fits.PrimaryHDU(frames[0])] + \
                       [fits.ImageHDU(frames[i]) \
                        for i in range(1, num_layers)])

    if filename is None:
        return hdu
    else:
        hdu.writeto(filename, overwrite=True, checksum=True)


def install_noise_cube(n=9):
    """
    Install a noise cube in the package directory

    Parameters
    ----------
    n : int, optional
        number of layers.

    Warning
    -------
    Each layer is ~64MB, default is 9 layers (~600MB). If you have less than
    1 GB on the drive where your Python installation is. Be careful!
    """

    if sys.version_info.major >= 3:
        print("WARNING - this process can take up to 10 minutes. Fear not!")
        hdu = make_noise_cube(n, filename=None)
        filename = find_file("FPA_noise.fits")
        hdu.writeto(filename, overwrite=True, checksum=True)
        print("Saved noise cube with", n, "layers to the package directory:")
        print(filename)
    else:
        print("Sorry, but this only works in Python 3 and above. \
               See the ScopeSim FAQs for work-around options")
