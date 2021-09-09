"""SpectralTraceList and SpectralTrace for the METIS LM spectrograph"""
import numpy as np

from astropy.io import ascii as ioascii   # ..todo: remove

from astropy.table import Table

from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace

class MetisLMSSpectralTraceList(SpectralTraceList):
    """
    SpectralTraceList for the METIS LM spectrograph
    """

    def __init__(self, **kwargs):
        # Set class parameters before initialisation
        # Call __init__ not with kwargs but with class parameters?
        pars = {
            "wavelen": kwargs['wavelen'],   # ..todo: move to yaml
            "naxis1": 122,
            "nslice": 28,
            "slicewidth": 0.0207,  # arcsec
            "pixscale": 0.0082,    # arcsec
            "grat_spacing": 18.2,
            "plate_scale": 0.303,
            }
        pars.update(kwargs)
        super().__init__(**pars)

        self.wavelen = kwargs['wavelen'] # ..todo: move to yaml
        self.nslice = 28
        # Parameters for METIS LMS layout
        # ..todo: Get this from layout file
        # ..todo: Update values
        self.naxis1 = 122
        self.nslice = 28
        self.slicewidth = 0.0207   # arcsec
        self.pixscale = 0.0082     # arcsec
        self.grat_spacing = 18.2
        self.plate_scale = 0.303

        # field of view of the instrument
        # ..todo: get this from aperture list
        self.view = np.array([self.naxis1 * self.pixscale,
                              self.nslice * self.slicewidth])


        # ..todo: fill with life
        #if self._file is not None:
        #    self.make_spectral_traces()

    def _class_params(self):
        """Parameters that are specific to the subclass"""
        params = {}
        self.meta.update(params)

    def make_spectral_traces(self):
        """
        Compute the transformations by interpolation
        """
        #nslice = len(self._file['Aperture List'].data)
        # determine echelle order and angle from specified wavelength
        tempres = self.angle_from_lambda()
        self.order = tempres['Ord']
        self.angle = tempres['Angle']

        spec_traces = dict()
        for sli in np.arange(self.meta['nslice']):
            slicename = "Slice " + str(sli + 1)
            spec_traces[slicename] = MetisLMSSpectralTrace(
                self._file["Polynomial Coefficients"],
                slice=sli, order=self.order, angle=self.angle,
                wavelen=self.meta['wavelen'])
            # ..todo: maybe pass self.meta rather than individual params?

        self.spectral_traces = spec_traces

    def angle_from_lambda(self):
        """
        Determine optimal echelle rotation angle for wavelength

        Parameters
        ----------
        lambda : float
                central wavelength in microns
        grat_spacing : float
                grating rule spacing in microns

        Returns
        -------
        a `dict` with entries
        - `Ord`: echelle order
        - `Angle`: grism angle
        - `Phase`: phase
        """
        lam = self.meta['wavelen']
        grat_spacing = self.meta['grat_spacing']

        # Read wcal extension of layout file
        # ..todo: switch to extension of FITS file
        wcal_file = "lms_dist_wcal.txt"
        wcal = ioascii.read(wcal_file, comment="^#", format="csv")

        # Compute angles, determine which order gives angle closest to zero
        angles = wcal['c0'] * lam + wcal['c1']
        imin = np.argmin(np.abs(angles))

        # Extract parameters
        order = wcal['Ord'][imin]
        angle = angles[imin]

        # Compute the phase corresponding to the wavelength
        phase = lam * order / (2 * grat_spacing)

        return {"Ord": order, "Angle": angle, "Phase": phase}


class MetisLMSSpectralTrace(SpectralTrace):
    """
    SpectralTrace for the METIS LM spectrograph
    """

    # ..todo:
    # In SpectralTrace, move all action to a method make_trace
    # Here, implement make_trace as well

    def __init__(self, polyhdu, **kwargs):
        super().__init__(polyhdu, **kwargs)
        self.meta['description'] = "Slice " + str(self.meta['slice'] + 1)


    def compute_interpolation_functions(self):
        self.matrices = self.get_matrices()


    def get_matrices(self):
        """
        Extract matrix from lms_dist_poly.txt

        Evaluate polynomial to obtain matrices A, B, AI and BI at grism angle
        given echelle order and slice number

        Parameters
        ----------
        order : int
            Echelle order
        spslice : int
            Slice number
        angle : float
            Grism angle in degrees

        Returns
        -------
        dict of four np.arrays of shape (4, 4) each
        """
        spslice = self.meta['slice']
        order = self.meta['order']
        angle = self.meta['angle']
        matnames = ['A', 'B', 'AI', 'BI']
        matrices = dict()

        poly = self.table
        for matid in range(4):
            select = ((poly['Ord'] == order) *
                      (poly['Sli'] == spslice) *
                      (poly['Mat'] == matid))
            if not np.any(select):
                raise KeyError("Combination of Order, Slice not found")

    # ..todo: use filename and instantiate the effect from the fits file
    #         Can Effect/DataContainer deal with multi-extension files?

    def __repr__(self):
        msg = '<MetisLMSSpectralTrace> "{}" : {} um : Order {} : Angle {}'\
            ''.format(self.meta["description"],
                      self.meta["wavelen"],
                      self.meta["order"],
                      self.meta["angle"])
        return msg
