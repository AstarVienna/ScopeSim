"""
NGHXRG by Bernard Rauscher  
see the paper: http://arxiv.org/abs/1509.06264  
downloaded from: http://jwst.nasa.gov/publications.html
"""

# dependencies include: astropy, numpy, scipy [, datetime, warnings, os]
import os
import datetime
# import warnings

import numpy as np
from scipy.ndimage.interpolation import zoom
from astropy.io import fits
from astropy.stats.funcs import median_absolute_deviation as mad

# import matplotlib.pyplot as plt # Handy for debugging

#warnings.filterwarnings('ignore')

__all__ = []

class HXRGNoise:
    """
    A class to generate HxRG noise frames
    
    HXRGNoise is a class for making realistic Teledyne HxRG system
    noise. The noise model includes correlated, uncorrelated,
    stationary, and non-stationary components. The default parameters
    make noise that resembles Channel 1 of JWST NIRSpec. NIRSpec uses
    H2RG detectors. They are read out using four video outputs at
    100.000 pix/s/output.
    """

    # These class variables are common to all HxRG detectors
    nghxrg_version = 2.3 # Software version

    def __init__(self, naxis1=None, naxis2=None, naxis3=None, n_out=None,
                 dt=None, nroh=None, nfoh=None, pca0_file=None, verbose=False,
                 reverse_scan_direction=False,
                 reference_pixel_border_width=None):
        """
        Simulate Teledyne HxRG+SIDECAR ASIC system noise.

        Parameters
        ----------
        naxis1 : int
            X-dimension of the FITS cube
        naxis2 : int
            Y-dimension of the FITS cube
        naxis3 : int
            Z-dimension of the FITS cube (number of up-the-ramp samples)
        n_out : int
            Number of detector outputs
        nfoh : int
            New frame overhead in rows. This allows for a short wait at the end 
            of a frame before starting the next one.
        nroh : int
            New row overhead in pixels. This allows for a short
            wait at the end of a row before starting the next one.
        dt : int
            Pixel dwell time in seconds
        pca0_file : str
            Name of a FITS file that contains PCA-zero
        verbose : bool
            Enable this to provide status reporting
        reference_pixel_border_width : int 
            Width of reference pixel border around image area
        reverse_scan_direction : bool
            Enable this to reverse the fast scanner readout directions. This
            capability was added to support Teledyne's programmable fast scan
            readout directions. The default setting =False corresponds to
            what HxRG detectors default to upon power up.
        """

        # ======================================================================
        #
        # DEFAULT CLOCKING PARAMETERS
        #
        # The following parameters define the default HxRG clocking pattern. The
        # parameters that define the default noise model are defined in the
        # mknoise() method.
        #
        # ======================================================================

        # Default clocking pattern is JWST NIRSpec
        self.naxis1    = 2048  if naxis1   is None else int(naxis1)
        self.naxis2    = 2048  if naxis2   is None else int(naxis2)
        self.naxis3    = 1     if naxis3   is None else int(naxis3)
        self.n_out     = 4     if n_out    is None else int(n_out)
        self.dt        = 1.e-5 if dt       is None else dt
        self.nroh      = 12    if nroh     is None else int(nroh)
        self.nfoh      = 1     if nfoh     is None else int(nfoh)
        self.reference_pixel_border_width = 4 \
                                            if reference_pixel_border_width is \
                                            None else reference_pixel_border_width

        # Initialize PCA-zero file and make sure that it exists and is a file
        self.pca0_file = os.getenv('NGHXRG_HOME')+'/nirspec_pca0.fits' if \
                         pca0_file is None else pca0_file
        if os.path.isfile(self.pca0_file) is False:
            print('There was an error finding pca0_file! Check to be')
            print('sure that the NGHXRG_HOME shell environment')
            print('variable is set correctly and that the')
            print('$NGHXRG_HOME/ directory contains the desired PCA0')
            print('file. The default is nirspec_pca0.fits.')
            os.sys.exit()


        # ======================================================================

        # Configure status reporting
        self.verbose = verbose

        # Configure readout direction
        self.reverse_scan_direction = reverse_scan_direction

        # Compute the number of pixels in the fast-scan direction per
        # output
        self.xsize = self.naxis1 // self.n_out

        # Compute the number of time steps per integration, per
        # output
        self.nstep = (self.xsize+self.nroh) * (self.naxis2+self.nfoh)\
                     * self.naxis3

        # For adding in ACN, it is handy to have masks of the even
        # and odd pixels on one output neglecting any gaps
        self.m_even = np.zeros((self.naxis3,self.naxis2,self.xsize))
        self.m_odd = np.zeros_like(self.m_even)
        for x in np.arange(0,self.xsize,2):
            self.m_even[:,:self.naxis2,x] = 1
            self.m_odd[:,:self.naxis2,x+1] = 1
        self.m_even = np.reshape(self.m_even, np.size(self.m_even))
        self.m_odd = np.reshape(self.m_odd, np.size(self.m_odd))

        # Also for adding in ACN, we need a mask that point to just
        # the real pixels in ordered vectors of just the even or odd
        # pixels
        self.m_short = np.zeros((self.naxis3, self.naxis2+self.nfoh, \
                                      (self.xsize+self.nroh)//2))
        self.m_short[:,:self.naxis2,:self.xsize//2] = 1
        self.m_short = np.reshape(self.m_short, np.size(self.m_short))

        # Define frequency arrays
        self.f1 = np.fft.rfftfreq(self.nstep)   # Frequencies for nstep elements
        self.f2 = np.fft.rfftfreq(2*self.nstep)  # ... for 2*nstep elements

        # Define pinkening filters. F1 and p_filter1 are used to
        # generate ACN. F2 and p_filter2 are used to generate 1/f noise.
        self.alpha = -1 # Hard code for 1/f noise until proven otherwise
        self.p_filter1 = np.sqrt(self.f1**self.alpha)
        self.p_filter2 = np.sqrt(self.f2**self.alpha)
        self.p_filter1[0] = 0.
        self.p_filter2[0] = 0.

        # Initialize pca0. This includes scaling to the correct size,
        # zero offsetting, and renormalization. We use robust statistics
        # because pca0 is real data
        hdu = fits.open(self.pca0_file)
        naxis1 = hdu[0].header['naxis1']
        naxis2 = hdu[0].header['naxis2']
        if (naxis1 != self.naxis1 or naxis2 != self.naxis2):
            zoom_factor = self.naxis1 / naxis1
            self.pca0 = zoom(hdu[0].data, zoom_factor, order=1, mode='wrap')
        else:
            self.pca0 = hdu[0].data
        self.pca0 -= np.median(self.pca0) # Zero offset
        self.pca0 /= (1.4826*mad(self.pca0)) # Renormalize

    def message(self, message_text):
        """
        Used for status reporting
        """
        if self.verbose is True:
            print('NG: ' + message_text + ' at DATETIME = ', \
                  datetime.datetime.now().time())

    def white_noise(self, nstep=None):
        """
        Generate white noise for an HxRG including all time steps
        (actual pixels and overheads).

        Parameters
        ----------
        nstep : int
            Length of vector returned
        """
        return(np.random.standard_normal(nstep))

    def pink_noise(self, mode):
        """
        Generate a vector of non-periodic pink noise.

        Parameters
        ----------
        mode : str
            Selected from {'pink', 'acn'}
        """

        # Configure depending on mode setting
        if mode is 'pink':
            nstep = 2*self.nstep
            f = self.f2
            p_filter = self.p_filter2
        else:
            nstep = self.nstep
            f = self.f1
            p_filter = self.p_filter1

        # Generate seed noise
        mynoise = self.white_noise(nstep)

        # Save the mean and standard deviation of the first
        # half. These are restored later. We do not subtract the mean
        # here. This happens when we multiply the FFT by the pinkening
        # filter which has no power at f=0.
        the_mean = np.mean(mynoise[:nstep//2])
        the_std = np.std(mynoise[:nstep//2])

        # Apply the pinkening filter.
        thefft = np.fft.rfft(mynoise)
        thefft = np.multiply(thefft, p_filter)
        result = np.fft.irfft(thefft)
        result = result[:nstep//2] # Keep 1st half

        # Restore the mean and standard deviation
        result *= the_std / np.std(result)
        result = result - np.mean(result) + the_mean

        # Done
        return(result)

    def mknoise(self, o_file, rd_noise=None, pedestal=None, c_pink=None,
                u_pink=None, acn=None, pca0_amp=None,
                reference_pixel_noise_ratio=None, ktc_noise=None,
                bias_offset=None, bias_amp=None):
        """
        Generate a FITS cube containing only noise.

        Parameters
        ----------
        o_file : str
            Output filename
        pedestal : float
            Magnitude of pedestal drift in electrons
        rd_noise : float
            Standard deviation of read noise in electrons
        c_pink : float
            Standard deviation of correlated pink noise in electrons
        u_pink : float
            Standard deviation of uncorrelated pink noise in electrons
        acn: float
            Standard deviation of alterating column noise in electrons
        pca0 : float
            Standard deviation of pca0 in electrons
        reference_pixel_noise_ratio : float
            Ratio of the standard deviation of
            the reference pixels to the regular
            pixels. Reference pixels are usually
            a little lower noise.
        ktc_noise : float
            kTC noise in electrons. Set this equal to
            sqrt(k*T*C_pixel)/q_e, where k is Boltzmann's
            constant, T is detector temperature, and C_pixel is
            pixel capacitance. For an H2RG, the pixel capacitance
            is typically about 40 fF.
        bias_offset : float
            On average, integrations stare here in electrons. Set this so that 
            all pixels are in range.
        bias_amp : float
            A multiplicative factor that we multiply PCA-zero by
            to simulate a bias pattern. This is completely
            independent from adding in "picture frame" noise.

        Notes
        -----
        Because of the noise correlations, there is no simple way to
        predict the noise of the simulated images. However, to a
        crude first approximation, these components add in
        quadrature.

        The units in the above are mostly "electrons". This follows convention
        in the astronomical community. From a physics perspective, holes are
        actually the physical entity that is collected in Teledyne's p-on-n
        (p-type implants in n-type bulk) HgCdTe architecture.

        """

        self.message('Starting mknoise()')

        # ======================================================================
        #
        # DEFAULT NOISE PARAMETERS
        #
        # These defaults create noise similar to that seen in the JWST NIRSpec.
        #
        # ======================================================================

        self.rd_noise  = 5.2      if rd_noise     is None else rd_noise
        self.pedestal  = 4        if pedestal     is None else pedestal
        self.c_pink    = 3        if c_pink       is None else c_pink
        self.u_pink    = 1        if u_pink       is None else u_pink
        self.acn       = .5       if acn          is None else acn
        self.pca0_amp  = .2       if pca0_amp     is None else pca0_amp

        # Change this only if you know that your detector is different from a
        # typical H2RG.
        self.reference_pixel_noise_ratio = 0.8 if \
            reference_pixel_noise_ratio is None else reference_pixel_noise_ratio

        # These are used only when generating cubes. They are
        # completely removed when the data are calibrated to
        # correlated double sampling or slope images. We include
        # them in here to make more realistic looking raw cubes.
        self.ktc_noise   = 29.   if ktc_noise   is None else ktc_noise
        self.bias_offset = 5000. if bias_offset is None else bias_offset
        self.bias_amp    = 500.  if bias_amp    is None else bias_amp

        # ======================================================================

        # Initialize the result cube. For up-the-ramp integrations,
        # we also add a bias pattern. Otherwise, we assume
        # that the aim was to simulate a two dimensional correlated
        # double sampling image or slope image.
        self.message('Initializing results cube')
        result = np.zeros((self.naxis3, self.naxis2, self.naxis1), \
                          dtype=np.float32)
        if self.naxis3 > 1:
            # Inject a bias pattern and kTC noise. If there are no reference pixels,
            # we know that we are dealing with a subarray. In this case, we do not
            # inject any bias pattern for now.
            if self.reference_pixel_border_width > 0:
                bias_pattern = self.pca0*self.bias_amp + self.bias_offset
            else:
                bias_pattern = self.bias_offset

            # Add in some kTC noise. Since this should always come out
            # in calibration, we do not attempt to model it in detail.
            bias_pattern += \
                         self.ktc_noise * \
                         np.random.standard_normal((self.naxis2, self.naxis1))

            # Ensure that there are no negative pixel values. Data cubes
            # are converted to unsigned integer before writing.
            bias_pattern = np.where(bias_pattern < 0, 0, bias_pattern)

            # Add in the bias pattern
            for z in np.arange(self.naxis3):
                result[z,:,:] += bias_pattern

        # Make white read noise. This is the same for all pixels.
        self.message('Generating rd_noise')
        w = self.reference_pixel_border_width # Easier to work with
        r = self.reference_pixel_noise_ratio  # Easier to work with
        for z in np.arange(self.naxis3):
            here = np.zeros((self.naxis2, self.naxis1))
            if w > 0: # Ref. pixel border exists
                # Add both reference and regular pixels
                here[:w, :] = r * self.rd_noise * \
                              np.random.standard_normal((w, self.naxis1))
                here[-w:, :] = r * self.rd_noise * \
                               np.random.standard_normal((w, self.naxis1))
                here[:, :w] = r * self.rd_noise * \
                              np.random.standard_normal((self.naxis1, w))
                here[:, -w:] = r * self.rd_noise * \
                               np.random.standard_normal((self.naxis1, w))
                # Make noisy regular pixels
                here[w:-w, w:-w] = self.rd_noise * \
                                  np.random.standard_normal((self.naxis2-2*w,
                                                             self.naxis1-2*w))
            else: # Ref. pixel border does not exist
                # Add only regular pixels
                here = self.rd_noise * np.random.standard_normal((self.naxis2,
                                                                  self.naxis1))
            # Add the noise in to the result
            result[z, :, :] += here

        # Add correlated pink noise.
        self.message('Adding c_pink noise')
        tt = self.c_pink * self.pink_noise('pink') # tt is a temp. variable
        tt = np.reshape(tt, (self.naxis3, self.naxis2+self.nfoh, \
                             self.xsize+self.nroh))[:,:self.naxis2,:self.xsize]
        for op in np.arange(self.n_out):
            x0 = op * self.xsize
            x1 = x0 + self.xsize
            if self.reverse_scan_direction is False:
                # Teledyne's default fast-scan directions
                if np.mod(op, 2) == 0:
                    result[:, :, x0:x1] += tt
                else:
                    result[:, :, x0:x1] += tt[:, :, ::-1]
            else:
                # Reverse the fast-scan directions.
                if np.mod(op,2)==1:
                    result[:,:,x0:x1] += tt
                else:
                    result[:,:,x0:x1] += tt[:,:,::-1]

        # Add uncorrelated pink noise. Because this pink noise is stationary and
        # different for each output, we don't need to flip it.
        self.message('Adding u_pink noise')
        for op in np.arange(self.n_out):
            x0 = op * self.xsize
            x1 = x0 + self.xsize
            tt = self.u_pink * self.pink_noise('pink')
            tt = np.reshape(tt, (self.naxis3, self.naxis2+self.nfoh, \
                            self.xsize+self.nroh))[:,:self.naxis2,:self.xsize]
            result[:, :, x0:x1] += tt

        # Add ACN
        self.message('Adding acn noise')
        for op in np.arange(self.n_out):

            # Generate new pink noise for each even and odd vector.
            # We give these the abstract names 'a' and 'b' so that we
            # can use a previously worked out formula to turn them
            # back into an image section.
            a = self.acn * self.pink_noise('acn')
            b = self.acn * self.pink_noise('acn')

            # Pick out just the real pixels (i.e. ignore the gaps)
            a = a[np.where(self.m_short == 1)]
            b = b[np.where(self.m_short == 1)]

            # Reformat into an image section. This uses the formula
            # mentioned above.
            acn_cube = np.reshape(np.transpose(np.vstack((a,b))),
                                  (self.naxis3,self.naxis2,self.xsize))

            # Add in the ACN. Because pink noise is stationary, we can
            # ignore the readout directions. There is no need to flip
            # acn_cube before adding it in.
            x0 = op * self.xsize
            x1 = x0 + self.xsize
            result[:,:,x0:x1] += acn_cube

        # Add PCA-zero. The PCA-zero template is modulated by 1/f.
        if self.pca0_amp > 0:
            self.message('Adding PCA-zero "picture frame" noise')
            gamma = self.pink_noise(mode='pink')
            zoom_factor = self.naxis2 * self.naxis3 / np.size(gamma)
            gamma = zoom(gamma, zoom_factor, order=1, mode='mirror')
            gamma = np.reshape(gamma, (self.naxis3, self.naxis2))
            for z in np.arange(self.naxis3):
                for y in np.arange(self.naxis2):
                    result[z, y, :] += self.pca0_amp*self.pca0[y, :]*gamma[z, y]

        # If the data cube has only 1 frame, reformat into a 2-dimensional
        # image.
        if self.naxis3 == 1:
            self.message('Reformatting cube into image')
            result = result[0, :, :]

        # If the data cube has more than one frame, convert to unsigned
        # integer
        if self.naxis3 > 1:
            self.message('Converting to 16-bit unsigned integer')
            result = result.astype('uint16')

        # Write the result to a FITS file
        self.message('Writing FITS file')
        hdu = fits.PrimaryHDU(result)
        hdu.header.append()
        hdu.header.append(('RD_NOISE', self.rd_noise, 'Read noise'))
        hdu.header.append(('PEDESTAL', self.pedestal, 'Pedestal drifts'))
        hdu.header.append(('C_PINK', self.c_pink, 'Correlated pink'))
        hdu.header.append(('U_PINK', self.u_pink, 'Uncorrelated pink'))
        hdu.header.append(('ACN', self.acn, 'Alternating column noise'))
        hdu.header.append(('PCA0', self.pca0_amp, \
                           'PCA zero, AKA picture frame'))
        #hdu.header['HISTORY'] = 'Created_by_NGHXRG_version_' \
        #                        + str(self.nghxrg_version)

        self.message('Exiting mknoise()')

        if o_file is not None:
            hdu.writeto(o_file, clobber='True')
        return result
