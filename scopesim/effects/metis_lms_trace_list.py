"""SpectralTraceList and SpectralTrace for the METIS LM spectrograph"""
import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table

from ..utils import from_currsys
from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace
from .spectral_trace_list_utils import Transform2D


class MetisLMSSpectralTraceList(SpectralTraceList):
    """
    SpectralTraceList for the METIS LM spectrograph
    """
    _class_params = {
        "naxis1": 122,
        "nslice": 28,
        "slicewidth": 0.0207,  # arcsec
        "pixscale": 0.0082,    # arcsec
        "grat_spacing": 18.2,
        "plate_scale": 0.303,
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        #self.params = {"wavelen": "!OBS.wavelen"}
        #self.params.update(kwargs)

        self.wavelen = self.meta['wavelen']

        # field of view of the instrument
        # ..todo: get this from aperture list
        self.view = np.array([self.meta['naxis1'] * self.meta['pixscale'],
                              self.meta['nslice'] * self.meta['slicewidth']])

        #if self._file is not None:
        #    print(self._file)
        #    self.make_spectral_traces()

    def make_spectral_traces(self):
        """
        Compute the transformations by interpolation
        """
        #nslice = len(self._file['Aperture List'].data)
        # determine echelle order and angle from specified wavelength
        tempres = self._angle_from_lambda()
        self.meta['order'] = tempres['Ord']
        self.meta['angle'] = tempres['Angle']

        spec_traces = dict()
        for sli in np.arange(self.meta['nslice']):
            slicename = "Slice " + str(sli + 1)
            spec_traces[slicename] = MetisLMSSpectralTrace(
                self._file,
                spslice=sli, params=self.meta)

        self.spectral_traces = spec_traces

    def _angle_from_lambda(self):
        """
        Determine optimal echelle rotation angle for wavelength
        """
        lam = from_currsys(self.meta['wavelen'])
        grat_spacing = self.meta['grat_spacing']
        wcal = self._file['WCAL'].data
        return echelle_setting(lam, grat_spacing, wcal)


class MetisLMSSpectralTrace(SpectralTrace):
    """
    SpectralTrace for the METIS LM spectrograph
    """
    _class_params = {
        "naxis1": 122,
        "nslice": 28,
        "slicewidth": 0.0207,  # arcsec
        "pixscale": 0.0082,    # arcsec
        "grat_spacing": 18.2,
        "plate_scale": 0.303,
    }

    def __init__(self, hdulist, spslice, params, **kwargs):
        polyhdu = hdulist['Polynomial coefficients']
        params.update(kwargs)
        params['aperture_id'] = spslice
        params['slice'] = spslice
        super().__init__(polyhdu, **params)
        self._file = hdulist
        self.meta['description'] = "Slice " + str(spslice + 1)
        self.meta.update(params)

    def fov_grid(self, fov_manager):
        """
        Provide information on the source space volume required by the effect

        Returns
        -------
        A dictionary with entries `wave_min` and `wave_max`, `x_min`, `y_min`,
        `x_max`, `y_max`. Spatial limits refer to the sky and are given in
        arcsec.
        """
        aperture = self._file['Aperture list'].data[self.meta['slice']]
        x_min = aperture['left']
        x_max = aperture['right']
        y_min = aperture['bottom']
        y_max = aperture['top']

        # ..todo: compute wave_min and wave_max
        return {"x_min": x_min, "x_max": x_max,
                "y_min": y_min, "y_max": y_max}


    def compute_interpolation_functions(self):
        """
        Define the transforms between (xi, lam) and (x, y).
        The LMS transforms actually operate on phase rather than
        wavelength, hence the necessity of defining pre- and
        posttransforms on the lam variable.
        """
        matrices = self.get_matrices()
        # matrices are transposed to align argument sequence
        # with the name of the functions
        self.xilam2x = Transform2D(matrices['A'].T,
                                   pretransform_x=self.sky2fp,
                                   pretransform_y=self.lam2phase)
        self.xilam2y = Transform2D(matrices['B'].T,
                                   pretransform_x=self.sky2fp,
                                   pretransform_y=self.lam2phase)
        self.xy2lam = Transform2D(matrices['AI'],
                                  posttransform=self.phase2lam)
        self.xy2xi = Transform2D(matrices['BI'],
                                 posttransform=self.fp2sky)


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

            subpoly = poly[select]
            thematrix = np.zeros((4, 4))
            for i in range(4):
                for j in range(4):
                    sel_ij = (subpoly['Row'] == i) * (subpoly['Col'] == j)
                    thematrix[i, j] = (subpoly['A11'][sel_ij] * angle**3 +
                                       subpoly['A12'][sel_ij] * angle**2 +
                                       subpoly['A21'][sel_ij] * angle +
                                       subpoly['A22'][sel_ij])
            matrices[matnames[matid]] = thematrix

        return matrices

    # ..todo: use filename and instantiate the effect from the fits file
    #         Can Effect/DataContainer deal with multi-extension files?

    def lam2phase(self, lam):
        """
        Convert wavelength to phase

        Phase is lam * order / (2 * grat_spacing).

        Parameters
        ----------
        lam : ndarray (float)
            wavelength (um)

        Returns
        -------
        Phase : ndarray
        """
        return self.meta['order'] * lam / (2 * self.meta['grat_spacing'])

    def phase2lam(self, phase):
        """
        Convert phase to wavelength

        Wavelength is phase * 2 * grat_spacing / order

        Parameters
        ----------
        phase : ndarray (float)
            phase (dimensionless)

        Returns
        -------
        wavelength : ndarray (um)
        """
        return 2 * self.meta['grat_spacing'] * phase / self.meta['order']

    def sky2fp(self, xi):
        """
        Convert position in arcsec to position in FP2
        """
        return xi / self.meta['plate_scale']

    def fp2sky(self, fp_x):
        """
        Convert position in FP2 to position on sky
        """
        return fp_x * self.meta['plate_scale']


    def __repr__(self):
        msg = '<MetisLMSSpectralTrace> "{}" : {} um : Order {} : Angle {}'\
            ''.format(self.meta["description"],
                      self.meta["wavelen"],
                      self.meta["order"],
                      self.meta["angle"])
        return msg


def echelle_setting(wavelength, grat_spacing, wcal_def):
    """
    Determine optimal echelle rotation angle for wavelength

    Parameters
    ----------
    lambda : float
            central wavelength in microns
    grat_spacing : float
            grating rule spacing in microns
    wcal_def: fits.TableHDU, fits.BinTableHDU, Table, str
            definition of the wavelength calibration parameters
            If str, interpreted as name of a fits file, with a
            table extension 'WCAL'.

    Returns
    -------
    a `dict` with entries
    - `Ord`: echelle order
    - `Angle`: grism angle
    - `Phase`: phase
    """
    if isinstance(wcal_def, fits.FITS_rec):
        wcal = wcal_def
    elif isinstance(wcal_def, (fits.TableHDU, fits.BinTableHDU)):
        # Read wcal extension of layout file
        wcal = wcal_def.data
    elif isinstance(wcal_def, Table):
        wcal = wcal_def
    elif isinstance(wcal_def, str):
        try:
            wcal = fits.getdata(wcal_def, extname='WCAL')
        except OSError:
            wcal = ioascii.read(wcal_def, comment="^#", format="csv")
    else:
        raise TypeError("wcal_def not in recognised format:", wcal_def)

    # Compute angles, determine which order gives angle closest to zero
    angles = wcal['c0'] * wavelength + wcal['c1']
    imin = np.argmin(np.abs(angles))

    # Extract parameters
    order = wcal['Ord'][imin]
    angle = angles[imin]

    # Compute the phase corresponding to the wavelength
    phase = wavelength * order / (2 * grat_spacing)

    return {"Ord": order, "Angle": angle, "Phase": phase}
