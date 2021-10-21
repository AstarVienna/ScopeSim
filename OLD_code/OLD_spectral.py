"""Classes for spectral curves"""
###############################################################################
# spectral
#
# DESCRIPTION
#
# "If you pay peanuts, you get monkeys"
#
# A spectral is the base class for either a transmission curve or an
# emission curve. The main attributes are 2 equal length arrays holding the
# centres of each wavelength bin and the corresponding value - an energy or a
# transmission factor [0-1]
#  - lam
#  - val
#
# spectral should be overloaded on the + and * operators. Although a rebin
# method would be good, this is unique to the subclasses. E.g. a
# Throughput curve rebin would involve averaging the "val" values, while a
# spectral rebin would involve summing up the "val" values.
#
# spectral also needs the file path of the data:
#  - filename
#
# A Throughput curve doesn't need anything else on top of the spectral,
# however the EmissionCurve must know which units are being used so that it
# can immediately convert the energy into photons. In order to do this the
# EmissionCurve needs the following extra info from the tests_commands dictionary
#  - spatial area [m2]
#  - angular area [arcsec2]
#  - integration time [s]
#
# As stated above each subclass should have its own rebin(lam_res) method
#
# Notes:
# All wavelength values are in [um]
# All other values are either transmission [0-1] or number of photons [>=0]
#
# Classes:
#  spectral(object)
#  - from_file(filename)
#  - from_list([ThroughputCurve])
#  - from_skycalc(filename)
#
# Subclasses:
# Emission(spectral)
# - rebin(lam)
# Throughput(spectral)
# - rebin(lam)
#
# Methods:
#
#
#

import os

from copy import deepcopy
import warnings

import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.io import ascii as ioascii  # 'ascii' redefines built-in
import astropy.table 
import yaml

from scopesim.utils import find_file

__all__ = ["TransmissionCurve", "EmissionCurve", "BlackbodyCurve", "UnityCurve",
           "get_sky_spectrum"]


class TransmissionCurve(object):
    """
    Very basic class to either read in a text file for a transmission curve
    or take two vectors to make a transmission curve

    Parameters
    ----------
    filename : str, optional
        The path to the file containing wavelength and transmission data
        where the first column is wavelength in [um] and the second is the
        transmission coefficient between [0,1].

        Alternatively this data can be passed directly. If filename is not
        provided, ``lam=`` and ``val=`` must be passed
    lam : array, optional
        [um] Wavelength bins in a 1D numpy array of length n
    val : array, optional
        [0 .. 1] The transmission coefficients

    Optional Parameters (**kwargs)
    ------------------------------
    lam_res : float
        [um] float with the desired spectral resolution. Default is 0.001.
    min_step : float
        [um] the minimum bin size used when resampling. Default is 1E-4
    wave_unit : str
        If the wavelength bins are in units other than micron. Default is "um"
    use_default_lam : bool
        Default is True. If True the curve is resmapled to a default range. This
        is useful for getting all curves on the same grid from the beginning
    default_lam : array
        Default is [0.3 .. 3.0] um with a step size of 0.001 um. A default
        wavelength range to avoid excessive resampling operations.
    airmass : float
        If transmission coefficients for more than one airmass are in ``filename``

    Returns
    -------
    TransmissionCurve object

    """
    def __init__(self, filename=None, lam=None, val=None, **kwargs):
        # default parameters. Can be overridden with kwargs.
        self.params = {"filename"    : filename,
                       "lam"         : lam,
                       "val"         : val,
                       "lam_res"     : 0.001,
                       "Type"        : "Transmission",
                       "min_step"    : 1E-4,
                       "wave_unit"    : u.um,
                       "use_default_lam" : True,
                       "on_default_lam" : False,
                       "default_lam" : np.arange(0.3, 3.0, 0.01),
                       "airmass"     : None}

        self.params.update(kwargs)

        self.info = dict([])
        self.info["Type"] = self.params["Type"]

        self.lam_orig, self.val_orig = self._get_data()
        self.lam_orig *= (1 * self.params["wave_unit"]).to(u.um).value

        self.lam = self.lam_orig
        self.val = self.val_orig

        # Resample to a default regular grid.
        # Not all of the filter curve .dat files have regular bin spacing and
        # so it is impossible to define a "lam_res" for those curves.
        # This is needed for the EmissionCurve method
        # "photons_in_range(lam_min, lam_max)"
        if self.params["Type"] == "Emission":
            self.resample(self.params["lam_res"], action="sum",
                          use_default_lam=self.params["use_default_lam"])
        else:
            self.resample(self.params["lam_res"], action="average",
                          use_default_lam=self.params["use_default_lam"])


    def __str__(self):
        return "Spectral curve:\n" + str(self.info)


    def __repr__(self):
        mask = [0, 1, 2], [-3, -2, -1]
        return self.info["Type"] + "Curve \n" + str(self.val[mask[0]])[:-1] + \
               " ..." + str(self.val[mask[1]])[1:]


    def _get_data(self):
        """
        Get the wavelength and value vectors from the input parameters

        Returns
        -------
        lam, val : 1D array

        """

        if self.params["lam"] is not None and self.params["val"] is not None:
            lam = self.params["lam"]
            val = self.params["val"]

        # test if it is a skycalc file
        elif self.params["filename"] is not None:
            filename = find_file(self.params["filename"])

            if ".fits" in filename:
                hdr = fits.getheader(filename)
                if any(["SKYCALC" in hdr[i] for i in range(len(hdr))
                        if isinstance(hdr[i], str)]):
                    if self.params["Type"] == "Emission":
                        lam = fits.getdata(filename)["lam"]
                        val = fits.getdata(filename)["flux"]
                    else:
                        lam = fits.getdata(filename)["lam"]
                        val = fits.getdata(filename)["trans"]
                else:
                    data = fits.getdata("../data/skytable.fits")
                    lam = data[data.columns[0].name]
                    val = data[data.columns[1].name]

            elif self.params["airmass"] is not None:
                lam, val = get_sky_spectrum(filename,
                                            airmass=self.params["airmass"])

            else:
                data = ioascii.read(filename)
                lam = data[data.colnames[0]]
                val = data[data.colnames[1]]
        else:
            raise ValueError("Please pass either filename or lam/val keywords")

        return lam, val


    def resample(self, bins, action="average", use_edges=False, min_step=None,
                 use_default_lam=False):
        """
        Resamples both the wavelength and value vectors to an even grid.
        In order to avoid losing spectral information, the TransmissionCurve
        resamples down to a resolution of 'min_step' (default: 0.01nm)
        before resampling again up to the given sampling vector defined by
        'bins'.

        Parameters
        ----------
        bins : float or array of floats
            [um]: float - taken to mean the width of bins on an even grid
                  array - the centres of the spectral bins
                        - the edges of the spectral bins if use_edges = True
        action : str, optional
            ['average','sum'] How to rebin the spectral curve. If 'sum',
            then the curve is normalised against the integrated value of
            the original curve. If 'average', the average value per bin
            becomes the value for each bin.
        use_edges : bool, optional
            [False, True] True if the array passed in 'bins' describes
            the edges of the wavelength bins. Default is False
        min_step : float, optional
            [um] default=1E-4, the step size for the down-sample
        use_default_lam : bool, optional
            Default is False. If True, ``bins`` is ignored and the default
            wavelength range is used as the resampling grid.
        """

        self.params["use_default_lam"] = use_default_lam
        if min_step is not None:
            self.params["min_step"] = min_step


        # No need to resample if the curve is already on the default grid
        if self.params["on_default_lam"] and self.params["use_default_lam"]:
            return

        #####################################################
        # Work out the irregular grid problem while summing #
        #####################################################

        tmp_x = np.arange(self.lam_orig[0], self.lam_orig[-1],
                          self.params["min_step"])
        tmp_y = np.interp(tmp_x, self.lam_orig, self.val_orig)

        #print("resampling", len(self.lam_orig))

        # use_default_lam overrides the bins argument
        if self.params["use_default_lam"]:
            lam_tmp = self.params["default_lam"]
            self.params["on_default_lam"] = True
        else:
            # if bins is a single number, use it as the bin width
            # else as the bin centres
            if not hasattr(bins, "__len__"):
                lam_tmp = np.arange(self.lam_orig[0],
                                    self.lam_orig[-1] + 1E-7, bins)
            else:
                lam_tmp = bins


        lam_res = lam_tmp[1] - lam_tmp[0]
        if self.params["min_step"] >= lam_res:
            warnings.warn("min_step > resample resolution. Can't resample")

        # define the edges and centres of each wavelength bin
        if use_edges:
            lam_bin_edges = lam_tmp
            lam_bin_centers = 0.5 * (lam_tmp[1:] + lam_tmp[:-1])
        else:
            lam_bin_edges = np.append(lam_tmp - 0.5*lam_res,
                                      lam_tmp[-1] + 0.5*lam_res)
            lam_bin_centers = lam_tmp

        # here is the assumption of a regular grid - see res_tmp
        val_tmp = np.zeros((len(lam_bin_centers)))

        for i in range(len(lam_bin_centers)):

            mask_i = np.where((tmp_x > lam_bin_edges[i]) *
                              (tmp_x < lam_bin_edges[i+1]))[0]

            if np.sum(mask_i) > 0 and action == "average":
                val_tmp[i] = np.average(tmp_y[mask_i[0]:mask_i[-1]])

            elif np.sum(mask_i) > 0 and action == "sum":
                # FIXED. THE SUMMING ISSUE. TEST IT         #
                # Tested - the errors are on the 0.1% level #
                val_tmp[i] = np.trapz(tmp_y[mask_i[0]:mask_i[-1]])
            else:
                val_tmp[i] = 0

        # The summing issue - assuming we want to integrate along the curve,
        # i.e. count all the photons in a new set of bins, we need to integrate
        # along the well-sampled (1E-5um) curve. However the above line of code
        # using np.interp doesn't take into account the new bin width when
        # resampling down to 1E-5um. I account for this summing up all the
        # photons in the original data set and normalising the new 1E-5 bin
        # data set to have the same amount.
        if action == "sum" and np.sum(val_tmp) != 0:
            val_tmp *= (np.sum(self.val_orig) / np.sum(val_tmp))

        self.lam = lam_tmp
        self.val = val_tmp
        self.res = lam_res
        self.params["lam_res"] = self.res


    def normalize(self, val=1., mode='integral'):
        """
        Normalize the spectral curve

        Parameters
        ----------
        val : float, optional
            The value to normalise to. Default is 1.
        mode : str, optional
            - "integral" normalizes the integral over the defined
               wavelength range to val (default: 1.)
            - "maximum" normalizes the maximum over the defined
               wavelength range to val (default: 1.)
        """
        if mode.lower() == 'integral':
            self.val = val * self.val / np.trapz(self.val, self.lam)
        elif mode.lower() == 'maximum':
            self.val = self.val / self.val.max() * val
        else:
            errorstr = "Unknown normalization mode: {0}. No action taken."
            raise ValueError(errorstr.format(mode))


    def plot(self, **kwargs):
        """
	Plot the transmission curve on the current axis

        The method accepts matplotlib.pyplot keywords.
        """
        import matplotlib.pyplot as plt
        plt.plot(self.lam, self.val, **kwargs)

    def filter_info(self):
        """
        Returns the filter properties as a dictionary

        Examples:
        ---------
        Creating a Table with fake values::

            >>> meta_dict =  {"comments": {"author : me" , "source : me", "date : today", "status : ready","center : 0.0",  "width : 1.1"}}
            >>> T = astropy.table.Table(data=[ [1.1,1.2,1.3], [0.1,0.2,0.3] ], names=("wavelength","transmission"), meta=meta_dict,copy=True)
            >>> T.write("tmp_table.dat",format="ascii",fast_writer=False)

        Reading the transmission curve:: 
            
            >>> Tc = TransmissionCurve("tmp_table.dat")
            >>> Tc.filter_info()

            {'author': 'me',
             'width': 1.1,
             'status': 'ready',
             'date': 'today',
             'source': 'me',
             'center': 0.0,
             'filename': 'tmp_table.dat'}

        Deleting the table::

            >>> os.remove('tmp_table.dat')

        """

        tbl = astropy.table.Table.read(self.params["filename"],format="ascii",header_start=-1)
        
        meta = tbl.meta      

        if not meta:
            cmts_dict = {"comments": ""}

        else:
            cmts_list = meta["comments"]
            cmts_str  = "\n".join(cmts_list)
            cmts_dict = yaml.full_load(cmts_str)
            if type(cmts_dict) is str:
                cmts_dict={"comments":cmts_dict}

        cmts_dict["filename"] = self.params["filename"]
        return cmts_dict


    def filter_table(self):
        """
        Returns the filter properties as a astropy.table

        Notes
        -----
        ONLY works if filter files have the ScopeSim header format

        The following keywords should be in the header::

            author
            source
            date_created
            date_modified
            status 
            type
            center
            width
            blue_cutoff
            red_cutoff

        """
        cmts_dict = self.filter_info()
        
        filter_table = astropy.table.Table()
        keys = [k for k in cmts_dict.keys()]

        req_keys = ['filename', 'center', 'width', 'blue_cutoff', 'red_cutoff', 
                    'author', 'source', 'date_created', 'date_modified', 'status', 'type']
	
        if np.all([k in req_keys for k in keys]):
            for keyword in req_keys:
                col = astropy.table.Column(name=keyword, data=(cmts_dict[keyword],))
                filter_table.add_column(col)
        else:
            raise ValueError(self.params["filename"] + " is not a ScopeSim filter")
        return filter_table


    def __len__(self):
        return len(self.val)


    def __getitem__(self, i):
        return self.val[i], self.lam[i]


    def __array__(self):
        return self.val


    def __pow__(self, n):
        """
        TransmissionCurve.val to the power of n.
        """
        tcnew = deepcopy(self)
        tcnew.val = tcnew.val ** n

        tcnew.lam_orig = tcnew.lam
        tcnew.val_orig = tcnew.val

        return tcnew


    def __mul__(self, tc):
        """
        Product of a TransmissionCurve with a scalar or another Curve
        If tc is a TransmissionCurve and does not have the same lam, it is
        resampled first.
        """

        tcnew = deepcopy(self)

        # EmissionCurve takes precedence over TransmissionCurve
        if not isinstance(self, EmissionCurve) and isinstance(tc, EmissionCurve):
            return tc * tcnew

        if not hasattr(tc, "val"):
            tcnew.val *= tc
        else:
            # Would np.allclose be better?
            if not np.all(self.lam == tc.lam):
                tc.resample(self.lam)
            tcnew.val *= tc.val

        tcnew.lam_orig = tcnew.lam
        tcnew.val_orig = tcnew.val

        return tcnew


    def __add__(self, tc):
        """
        Addition to a TransmissionCurve of a scalar or another Curve.
        If tc is a TransmissionCurve and does not have the same lam, it is
        resampled first.
        """
        tcnew = deepcopy(self)

        if not hasattr(tc, "val"):
            tcnew.val += tc
        else:
            ### Would np.allclose be better?
            if not np.all(self.lam == tc.lam):
                tc.resample(self.lam)
            tcnew.val += tc.val

        tcnew.lam_orig = tcnew.lam
        tcnew.val_orig = tcnew.val

        return tcnew


    def __sub__(self, tc):
        return  self.__add__(-1 * tc)


    def __rmul__(self, x):
        return self.__mul__(x)


    def __radd__(self, x):
        return self.__add__(x)


    def __rsub__(self, x):
        self.__mul__(-1)
        return self.__add__(x)


    def __imul__(self, x):
        return self.__mul__(x)


    def __iadd__(self, x):
        return self.__add__(x)


    def __isub__(self, x):
        return self.__sub__(x)



class EmissionCurve(TransmissionCurve):
    """
    Class for emission curves

    Create an emission curve from a file or from wavelength and flux vectors.
    In the latter case, the keywords `lam` and `val` have to be specified.

    Parameters
    ==========
    - filename: string with the path to the transmission curve file where
                the first column is wavelength in [um] and the second is the
                transmission coefficient between [0,1]
    - lam: [um] 1D numpy array of length n
    - val: 1D numpy array of length n
    - res: [um] float with the desired spectral resolution
    - pix_res: [arcsec] float of int for the field of view for each pixel
    - area: [m2] float or int for the collecting area of M1
    - units: string or astropy.unit for calculating the number of photons
             per voxel

    Return values are in [ph/s/voxel]

    Examples
    --------
    ::

        >>> from scopesim.spectral import EmissionCurve
        >>>
        >>> ec_1 = EmissionCurve("emission_curve.dat")
        >>> lam = np.arange(0.7, 1.5, 0.05)
        >>> flux = 1. - 0.2 * wave**2   # power-law spectrum
        >>> ec_2 = EmissionCurve(lam=lam, val=flux)
    """
    def __init__(self, filename=None, **kwargs):
        default_params = {"exptime" :1,
                          "pix_res" :0.004,
                          "area"    :978,
                          "units"   :"ph/(s m2 um arcsec2)"}
        if "units" not in kwargs.keys():
            warnings.warn("""
            No 'units' specified in EmissionCurve.
            Assuming ph/(s m2 micron arcsec2)""", RuntimeWarning)
        default_params.update(kwargs)

        if default_params["area"] < 1E-6:
            default_params["area"] = 1E-6

        if filename is not None:
            default_params["filename"] = filename

        super(EmissionCurve, self).__init__(Type="Emission", **default_params)
        self.factor = 1
        self.convert_to_photons()


    def resample(self, bins, action="average", use_edges=False, min_step=None,
                 use_default_lam=False):
        """Rebin an emission curve"""
        super(EmissionCurve, self).resample(bins=bins, action=action,
                                            use_edges=use_edges, min_step=min_step,
                                            use_default_lam=use_default_lam)


    def convert_to_photons(self):
        """Do the conversion to photons/s/voxel by using the val_unit, lam, area
        and exptime keywords. If not given, make some assumptions.
        """
        self.params["units"] = u.Unit(self.params["units"])
        bases = self.params["units"].bases

        factor = 1. * self.params["units"]

        # The delivered EmissionCurve should be in ph/s/voxel
        if u.s  not in bases:
            factor /= self.params["exptime"] * u.s
        if u.m      in bases:
            factor *= self.params["area"] * u.m**2
        if u.arcsec in bases:
            factor *= (self.params["pix_res"] *u.arcsec)**2
        if u.micron in bases:
            factor *= self.params["lam_res"] * u.um

        self.params["units"] = factor.unit

        self.val *= factor.value
        self.factor = factor


    def photons_in_range(self, lam_min=None, lam_max=None):
        """
        Sum up the photons in the wavelength range [lam_min, lam_max]

        Parameters
        ==========
        - lam_min, lam_max: the wavelength limits

        Return values are in [ph/s/pixel]
        """
        if lam_min is None:
            lam_min = self.lam[0]
        if lam_max is None:
            lam_max = self.lam[-1]

        if lam_min > self.lam[-1] or lam_max < self.lam[0]:
            warnings.warn("wavelength limits outside wavelength range")
            photons = 0
        else:
            if (lam_max - lam_min) < 2 * self.res:
                zoom = 10
                lam_zoom = np.linspace(lam_min, lam_max, zoom)
                spec_zoom = np.interp(lam_zoom, self.lam, self.val) / zoom
                photons = np.sum(spec_zoom)
            else:
                mask = (self.lam >= lam_min) * (self.lam < lam_max)
                photons = np.sum(self.val[mask])

        return photons


class BlackbodyCurve(EmissionCurve):
    """
    Blackbody emission curve

    Parameters
    ----------
    lam : 1D np.ndarray
        [um] the centres of the wavelength bins
    temp : float
        [deg C] float for the average temperature of the blackbody
    pix_res : float, optional
        [arcsec] Default is 0.004. Field of view for each pixel
    area : float, optional
        [m2] Default is 978m2. Area of the emitting surface

    Returns
    -------
    EmissionCurve object with units of [ph/s/voxel], i.e. photons per second
        per wavelength bin per full area per pixel field of view
    """

    def __init__(self, lam, temp, **kwargs):
        self.params = {"pix_res": 0.004, "area": 978, "temperature": temp}
        self.params.update(kwargs)

        # FOr some wierd reason it returns NaNs for temp < 266
        # if self.params["area"] < 1E-6:
        #    self.params["area"] = 1E-6

        temp += 273.15

        lam_res = lam[1] - lam[0]
        edges = np.append(lam - 0.5 * lam_res, lam[-1] + 0.5 * lam_res)
        lam_res = edges[1:] - edges[:-1]

        # I is in W sr-1 m-3 : [erg/s/sr/cm2/um]
        exparg = c.h * c.c / (c.k_B * (temp * u.K) * (lam * u.um))

        if np.any(exparg.si > 500):  # Wien approximation to avoid overflow
            I = 2. * c.h * c.c**2 /(lam * u.um)**5 * np.exp(-exparg)
        else:                        # Full Planck formula
            I = 2. * c.h * c.c**2 / (lam * u.um)**5 / (np.exp(exparg) - 1.)

        I = I / u.sr    # make it explicit that this is per steradian

        # E is in W
        E = I * (self.params["area"] * u.m**2) * (lam_res * u.um) * \
                (self.params["pix_res"] * u.arcsec)**2
        # ph is in 1/s
        ph = E / (c.h * c.c / (lam * u.um))
        val = ph.si.value

        super(BlackbodyCurve, self).__init__(lam=lam, val=val, units="1/s",
                                             **self.params)
        self.info["Type"] = "Blackbody"


class UnityCurve(TransmissionCurve):
    """Constant transmission curve

    Parameters
    ==========
    - lam [um]: wavelength array
    - val: constant value of transmission (default: 1)
    """

    def __init__(self, lam=np.asarray([0.3, 3.0]), val=1, **kwargs):
        val = np.asarray([val]*len(lam))
        super(UnityCurve, self).__init__(lam=lam, val=val, **kwargs)



def get_sky_spectrum(fname, airmass, return_type=None, **kwargs):
    """
    Return a spectral curve for the sky for a certain airmass

    Parameters
    ----------
    fname : str
        the file containing the spectral curves
    airmass : float, optinal
        Default is 1.0. Acceptable values are between 1.0 and 3.0
    return_type : str, optional
        ["transmission", "emission", None] Default is None. A TransmissionCurve
        or EmissionCurve object will be returned if desired. If None two arrays
        are returned: (lam, val)
    **kwargs : optional
        kwargs are passed directly onto the TransmissionCurve or EmissionCurve
        classes

    Returns
    -------
    TransmissionCurve or EmissionCurve or (lam, val)
        By default lam is in [um] and val [ph/s/m2/um/arcsec2] if val is an
        emission spectrum

    Notes
    -----
    This function is designed to work with a table of values produced by SkyCalc
    The column names must begin with lambda and then columns must be named
    according to the airmass following this pattern: "X1.5" for an airmass of 1.5
    """

    if not os.path.exists(fname):
        raise OSError("File doesn't exist: " + fname)

    data = ioascii.read(fname)
    tbl_airmass = np.array([float(i[1:]) for i in data.colnames[2:]])

    lam = data[data.colnames[0]]

    if "X" + str(float(airmass)) in data.colnames:
        val = data["X" + str(float(airmass))]


    elif airmass > 1 and airmass < 3:
        i = np.where(tbl_airmass - airmass >= 0)[0][0]

        # get the nearest two columns to the given airmass
        x0, x1 = tbl_airmass[(i-1):(i+1)]
        w = np.sum(tbl_airmass[(i-1):(i+1)] * np.array([-1, 1]))

        # get the weights for summing the two columns
        f1, f0 = (airmass - x0, x1 - airmass) / w    # backwards for a reason!
        val = f0 * data["X" + str(x0)] + f1 * data["X" + str(x1)]
    else:
        print("Column not found")

    if return_type is not None:
        if "trans" in return_type.lower():
            return TransmissionCurve(lam=lam, val=val, **kwargs)
        elif "emis" in return_type.lower():
            return EmissionCurve(lam=lam, val=val, **kwargs)
    else:
        return lam, val
