# -*- coding: utf-8 -*-
"""
Adapted from: NGHXRG by Bernard Rauscher.

see the paper: http://arxiv.org/abs/1509.06264
downloaded from: http://jwst.nasa.gov/publications.html
"""

# pylint: disable=no-member

from pathlib import Path
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import zoom
from astropy.io import fits
from astropy.stats.funcs import median_absolute_deviation as mad

from ..utils import get_logger


logger = get_logger(__name__)


@dataclass
class HXRGNoise:
    """
    Simulate Teledyne HxRG+SIDECAR ASIC system noise.

    HXRGNoise is a class for making realistic Teledyne HxRG system noise.
    The noise model includes correlated, uncorrelated, stationary, and
    non-stationary components. The default parameters make noise that resembles
    Channel 1 of JWST NIRSpec. NIRSpec uses H2RG detectors. They are read out
    using four video outputs at 100.000 pix/s/output.

    The following parameters define the default HxRG clocking pattern.
    The parameters that define the default noise model are defined in the
    mknoise() method. Default clocking pattern is JWST NIRSpec.

    Parameters
    ----------
    shape : tuple[int, int, int]
        Shape of the FITS cube [Z, Y, X] (Z...number of up-the-ramp samples).
    n_out : int
        Number of detector outputs
    n_foh : int
        New frame overhead in rows. This allows for a short wait at the end
        of a frame before starting the next one.
    n_roh : int
        New row overhead in pixels. This allows for a short
        wait at the end of a row before starting the next one.
    pca0_file : str
        Name of a FITS file that contains PCA-zero
    reference_pixel_border_width : int
        Width of reference pixel border around image area
    reverse_scan_direction : bool
        Enable this to reverse the fast scanner readout directions. This
        capability was added to support Teledyne's programmable fast scan
        readout directions. The default setting =False corresponds to
        what HxRG detectors default to upon power up.
    """

    shape: tuple[int, int, int] = (1, 2048, 2048)
    n_out: int = 4
    n_roh: int = 12
    n_foh: int = 1
    reference_pixel_border_width: int = 4
    pca0_file: str = None
    reverse_scan_direction: bool = False

    def __post_init__(self):
        """Run after dataclass init."""
        # Initialize PCA-zero file and make sure that it exists and is a file
        if self.pca0_file is None:
            self.pca0_file = Path("./nirspec_pca0.fits")
        else:
            self.pca0_file = Path(self.pca0_file)
        # if not pca0_file.exists():
        #     print("There was an error finding pca0_file! Check to be")
        #     print("sure that the NGHXRG_HOME shell environment")
        #     print("variable is set correctly and that the")
        #     print("$NGHXRG_HOME/ directory contains the desired PCA0")
        #     print("file. The default is nirspec_pca0.fits.")

        # Configure readout direction
        self.reverse_scan_direction = self.reverse_scan_direction

        # Compute the number of pixels in the fast-scan direction per
        # output
        # self.xsize = self.shape[2] // self.n_out
        self.scan_shape = np.array([
            *self.shape[:2],
            self.shape[2] // self.n_out
        ])

        # Compute the number of time steps per integration, per
        # output
        self.nstep = ((self.xsize + self.n_roh) *
                      (self.shape[1] + self.n_foh) *
                      self.shape[0])

        # For adding in ACN, it is handy to have masks of the even
        # and odd pixels on one output neglecting any gaps
        m_even = np.zeros(self.scan_shape)
        m_odd = np.zeros_like(m_even)
        for x in np.arange(0, self.xsize, 2):
            m_even[:, :self.shape[1], x] = 1
            m_odd[:, :self.shape[1], x + 1] = 1
        m_even = np.reshape(m_even, np.size(m_even))
        m_odd = np.reshape(m_odd, np.size(m_odd))

        # Also for adding in ACN, we need a mask that point to just
        # the real pixels in ordered vectors of just the even or odd
        # pixels
        self.m_short = np.zeros((self.shape[0], self.shape[1] + self.n_foh,
                                 (self.xsize + self.n_roh) // 2))
        self.m_short[:, :self.shape[1], :self.xsize // 2] = 1
        self.m_short = np.reshape(self.m_short, np.size(self.m_short))

        # Initialize pca0. This includes scaling to the correct size,
        # zero offsetting, and renormalization. We use robust statistics
        # because pca0 is real data
        with fits.open(self.pca0_file) as hdul:
            naxis1 = hdul[0].header["naxis1"]
            naxis2 = hdul[0].header["naxis2"]
            if (naxis1 != self.shape[2] or naxis2 != self.shape[1]):
                zoom_factor = self.shape[2] / naxis1
                self.pca0 = zoom(hdul[0].data, zoom_factor,
                                 order=1, mode="wrap")
            else:
                self.pca0 = hdul[0].data

        self.pca0 -= np.median(self.pca0)  # Zero offset
        self.pca0 /= (1.4826 * mad(self.pca0))  # Renormalize

    def white_noise(self, nstep: int):
        """
        Generate white noise for an HxRG including all time steps.

        (actual pixels and overheads).

        Parameters
        ----------
        nstep : int
            Length of vector returned
        """
        return np.random.standard_normal(nstep)

    def pink_noise(self, mode: str):
        """
        Generate a vector of non-periodic pink noise.

        Parameters
        ----------
        mode : str
            Selected from {"pink", "acn"}
        """
        # Configure depending on mode setting

        # Define pinkening filters. F1 and p_filter1 are used to
        # generate ACN. F2 and p_filter2 are used to generate 1/f noise.
        alpha = -1  # Hard code for 1/f noise until proven otherwise

        if mode == "pink":
            nstep = 2 * self.nstep
            # Frequencies for 2*nstep elements
            freq = np.fft.rfftfreq(2 * self.nstep)
            p_filter = np.sqrt(freq**alpha)
        else:
            nstep = self.nstep
            # Frequencies for nstep elements
            freq = np.fft.rfftfreq(self.nstep)
            p_filter = np.sqrt(freq**alpha)
        p_filter[0] = 0.

        # Generate seed noise
        mynoise = self.white_noise(nstep)

        # Save the mean and standard deviation of the first half. These are
        # restored later. We don't subtract the mean here. This happens when we
        # multiply the FFT by the pinkening filter which has no power at f=0.
        the_mean = np.mean(mynoise[:nstep // 2])
        the_std = np.std(mynoise[:nstep // 2])

        # Apply the pinkening filter.
        thefft = np.fft.rfft(mynoise)
        thefft = np.multiply(thefft, p_filter)
        result = np.fft.irfft(thefft)
        result = result[:nstep // 2]  # Keep 1st half

        # Restore the mean and standard deviation
        result *= the_std / result.std()
        result -= result.mean() + the_mean

        return result

    def _x_slice(self, op):
        x0 = op * self.scan_shape[2]
        x1 = x0 + self.scan_shape[2]
        return slice(x0, x1)

    def _mk_rd_noise(self, result, rd_noise, ref_px_noise_ratio):
        """Make white read noise. This is the same for all pixels."""
        logger.debug("Generating rd_noise")
        w = self.reference_pixel_border_width  # Easier to work with
        r = ref_px_noise_ratio   # Easier to work with
        for z in np.arange(self.shape[0]):
            here = np.zeros(self.shape[1:])
            if w > 0:  # Ref. pixel border exists
                # Add both reference and regular pixels
                here[:w, :] = (r * rd_noise *
                               np.random.standard_normal((w, self.shape[2])))
                here[-w:, :] = (r * rd_noise *
                                np.random.standard_normal((w, self.shape[2])))
                here[:, :w] = (r * rd_noise *
                               np.random.standard_normal((self.shape[2], w)))
                here[:, -w:] = (r * rd_noise *
                                np.random.standard_normal((self.shape[2], w)))
                # Make noisy regular pixels
                here[w:-w, w:-w] = (rd_noise *
                                    np.random.standard_normal(
                                        (self.shape[1] - 2 * w,
                                         self.shape[2] - 2 * w))
                                    )
            else:  # Ref. pixel border does not exist
                # Add only regular pixels
                here = rd_noise * np.random.standard_normal(self.shape[1:])
            # Add the noise in to the result
            result[z] += here
        return result

    def _mk_c_pink_noise(self, result, c_pink):
        """Add correlated pink noise."""
        logger.debug("Adding c_pink noise")
        tt = c_pink * self.pink_noise("pink")
        tt = np.reshape(tt, self.scan_shape + (0, self.n_foh, self.n_roh)
                        )[:, :self.shape[1], :self.xsize]
        modnum = int(self.reverse_scan_direction)
        for op in np.arange(self.n_out):
            x_slice = self._x_slice(op)
            if op % 2 == modnum:
                result[..., x_slice] += tt
            else:
                result[..., x_slice] += tt[..., ::-1]
        return result

    def _mk_u_pink_noise(self, result, u_pink):
        """Add uncorrelated pink noise.

        Because this pink noise is stationary and different for each output,
        we don't need to flip it.
        """
        logger.debug("Adding u_pink noise")
        for op in np.arange(self.n_out):
            tt = u_pink * self.pink_noise("pink")
            tt = np.reshape(tt, self.scan_shape + (0, self.n_foh, self.n_roh)
                            )[:, :self.shape[1], :self.xsize]
            result[..., self._x_slice(op)] += tt
        return result

    def _mk_acn_noise(self, result, acn):
        """Add ACN."""
        logger.debug("Adding acn noise")
        for op in np.arange(self.n_out):
            # Generate new pink noise for each even and odd vector.
            # We give these the abstract names "a" and "b" so that we
            # can use a previously worked out formula to turn them
            # back into an image section.
            a = acn * self.pink_noise("acn")
            b = acn * self.pink_noise("acn")

            # Pick out just the real pixels (i.e. ignore the gaps)
            a = a[np.where(self.m_short == 1)]
            b = b[np.where(self.m_short == 1)]

            # Reformat into an image section. This uses the formula
            # mentioned above.
            acn_cube = np.reshape(np.transpose(np.vstack((a, b))), self.scan_shape)

            # Add in the ACN. Because pink noise is stationary, we can
            # ignore the readout directions. There is no need to flip
            # acn_cube before adding it in.
            result[..., self._x_slice(op)] += acn_cube
        return result

    def _mk_pac0(self, result, pca0_amp):
        """Add PCA-zero. The PCA-zero template is modulated by 1/f."""
        if pca0_amp > 0:
            logger.debug("Adding PCA-zero 'picture frame' noise")
            gamma = self.pink_noise(mode="pink")
            zoom_factor = self.shape[1] * self.shape[0] / np.size(gamma)
            gamma = zoom(gamma, zoom_factor, order=1, mode="mirror")
            gamma = np.reshape(gamma, self.shape[:2])

            # TODO: vectorize if possible...
            for z in range(self.shape[0]):
                for y in range(self.shape[1]):
                    result[z, y] += pca0_amp * self.pca0[y] * gamma[z, y]
        return result

    def mknoise(
            self,
            o_file: str,
            rd_noise: float = 5.2,
            c_pink: float = 3.,
            u_pink: float = 1.,
            acn: float = .5,
            pca0_amp: float = .2,
            ref_px_noise_ratio: float = 0.8,
            ktc_noise: float = 29.,
            bias_offset: float = 5000.,
            bias_amp: float = 500.,
    ):
        """
        Generate a FITS cube containing only noise.

        Parameters
        ----------
        o_file : str
            Output filename
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
        ref_px_noise_ratio : float
            Ratio of the standard deviation of the reference pixels to the
            regular pixels. Reference pixels are usually a little lower noise.
            Change this only if you know that your detector is different from a
            typical H2RG.

        ktc_noise : float
            kTC noise in electrons. Set this equal to sqrt(k*T*C_pixel)/q_e,
            where k is Boltzmann's constant, T is detector temperature, and
            C_pixel is pixel capacitance. For an H2RG, the pixel capacitance
            is typically about 40 fF.
        bias_offset : float
            On average, integrations stare here in electrons. Set this so that
            all pixels are in range.
        bias_amp : float
            A multiplicative factor that we multiply PCA-zero by to simulate a
            bias pattern. This is completely independent from adding in
            "picture frame" noise.

        Notes
        -----
        Because of the noise correlations, there is no simple way to predict
        the noise of the simulated images. However, to a crude first
        approximation, these components add in quadrature.

        The units in the above are mostly "electrons". This follows convention
        in the astronomical community. From a physics perspective, holes are
        actually the physical entity that is collected in Teledyne's p-on-n
        (p-type implants in n-type bulk) HgCdTe architecture.

        The parameters `ktc_noise`, `bias_offset` and `bias_amp` are used only
        when generating cubes. They are completely removed when the data are
        calibrated to correlated double sampling or slope images. We include
        them in here to make more realistic looking raw cubes.

        """
        logger.debug("Starting mknoise()")

        # Initialize the result cube. For up-the-ramp integrations,
        # we also add a bias pattern. Otherwise, we assume
        # that the aim was to simulate a two dimensional correlated
        # double sampling image or slope image.
        logger.debug("Initializing results cube")
        result = np.zeros(self.shape, dtype=np.float32)
        if self.shape[0] > 1:
            # Inject a bias pattern and kTC noise. If there are no reference
            # pixels, we know that we are dealing with a subarray.
            # In this case, we do not inject any bias pattern for now.
            if self.reference_pixel_border_width > 0:
                bias_pattern = self.pca0 * bias_amp + bias_offset
            else:
                bias_pattern = bias_offset

            # Add in some kTC noise. Since this should always come out
            # in calibration, we do not attempt to model it in detail.
            bias_pattern += (
                ktc_noise *
                np.random.standard_normal(self.shape[1:])
            )

            # Ensure that there are no negative pixel values. Data cubes
            # are converted to unsigned integer before writing.
            bias_pattern = np.where(bias_pattern < 0, 0, bias_pattern)

            # Add in the bias pattern
            for z in range(self.shape[0]):
                result[z] += bias_pattern

        result = self._mk_rd_noise(result, rd_noise, ref_px_noise_ratio)
        result = self._mk_c_pink_noise(result, c_pink)
        result = self._mk_u_pink_noise(result, u_pink)
        result = self._mk_acn_noise(result, acn)
        result = self._mk_pac0(result, pca0_amp)

        # If the data cube has only 1 frame, reformat into a 2-dimensional
        # image.
        if self.shape[0] == 1:
            logger.debug("Reformatting cube into image")
            result = result[0]

        # If the data cube has more than one frame, convert to unsigned
        # integer
        if self.shape[0] > 1:
            logger.debug("Converting to 16-bit unsigned integer")
            result = result.astype("uint16")

        # Write the result to a FITS file
        logger.debug("Writing FITS file")
        hdu = fits.PrimaryHDU(result)
        hdu.header.append()
        hdu.header.append(("RD_NOISE", rd_noise, "Read noise"))
        hdu.header.append(("C_PINK", c_pink, "Correlated pink"))
        hdu.header.append(("U_PINK", u_pink, "Uncorrelated pink"))
        hdu.header.append(("ACN", acn, "Alternating column noise"))
        hdu.header.append(("PCA0", pca0_amp,
                           "PCA zero, AKA picture frame"))

        logger.debug("Exiting mknoise()")

        if o_file is not None:
            hdu.writeto(o_file)
        return result
