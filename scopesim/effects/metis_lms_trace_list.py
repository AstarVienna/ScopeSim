# -*- coding: utf-8 -*-
"""SpectralTraceList and SpectralTrace for the METIS LM spectrograph."""

import copy
import warnings

import numpy as np
from scipy.interpolate import RectBivariateSpline

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u

from ..utils import from_currsys, find_file, quantify, get_logger
from .spectral_trace_list import SpectralTraceList
from .spectral_trace_list_utils import SpectralTrace
from .spectral_trace_list_utils import Transform2D
from .spectral_trace_list_utils import make_image_interpolations
from .apertures import ApertureMask
from .ter_curves import TERCurve
from ..base_classes import FieldOfViewBase, FOVSetupBase
from ..optics.fov import FieldOfView


logger = get_logger(__name__)


class MetisLMSSpectralTraceList(SpectralTraceList):
    """SpectralTraceList for the METIS LM spectrograph."""

    _class_params = {
        "naxis1": 122,
        "nslice": 28,
        "slicewidth": 0.0207,  # arcsec
        "pixscale": 0.0082,    # arcsec
        "grat_spacing": 18.2,
        "fp2_platescale": 0.303,
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # self.params = {"wavelen": "!OBS.wavelen"}
        # self.params.update(kwargs)

        self.wavelen = self.meta["wavelen"]

        # field of view of the instrument
        self.slicelist = self._file["Aperture List"].data

        # self.view = np.array([self.meta["naxis1"] * self.meta["pixscale"],
        #                       self.meta["nslice"] * self.meta["slicewidth"]])

        self.view = np.array([self.slicelist["right"].max() -
                              self.slicelist["left"].min(),
                              self.slicelist["top"].max() -
                              self.slicelist["bottom"].min()])

        # for sli, spt in enumerate(self.spectral_traces.values()):
        #     spt.meta["xmin"] = self.slicelist["left"][sli]
        #     spt.meta["xmax"] = self.slicelist["right"][sli]
        #     spt.meta["ymin"] = self.slicelist["bottom"][sli]
        #     spt.meta["ymax"] = self.slicelist["top"][sli]

        # if self._file is not None:
        #     self.make_spectral_traces()

    def apply_to(self, obj, **kwargs):
        """See parent docstring."""
        if isinstance(obj, FOVSetupBase):
            # Create a single volume that covers the aperture and
            # the maximum wavelength range of LMS
            volumes = [spectral_trace.fov_grid()
                       for spectral_trace in self.spectral_traces.values()]

            wave_min = min(vol["wave_min"] for vol in volumes)
            wave_max = max(vol["wave_max"] for vol in volumes)
            extracted_vols = obj.extract(axes=["wave"],
                                         edges=([[wave_min, wave_max]]))
            obj.volumes = extracted_vols

        if isinstance(obj, FieldOfViewBase):
            # Application to field of view
            if obj.hdu is not None and obj.hdu.header["NAXIS"] == 3:
                obj.cube = obj.hdu
            elif obj.hdu is None and obj.cube is None:
                obj.cube = obj.make_cube_hdu()

            fovcube = obj.cube.data
            n_z, n_y, n_x = fovcube.shape
            fovwcs = WCS(obj.cube.header)
            # Make this linear to avoid jump at RA 0 deg
            fovwcs.wcs.ctype = ["LINEAR", "LINEAR", fovwcs.wcs.ctype[2]]
            fovwcs_spat = fovwcs.sub(2)
            ny_slice = self.meta["slice_samples"]

            # Spatial pixel coordinates
            xslice, yslice = np.meshgrid(np.arange(n_x),
                                         np.arange(ny_slice))

            fovimage = np.zeros((obj.detector_header["NAXIS2"],
                                 obj.detector_header["NAXIS1"]),
                                dtype=np.float32)

            for sptid, spt in self.spectral_traces.items():
                ymin = spt.meta["fov"]["y_min"]
                ymax = spt.meta["fov"]["y_max"]

                slicewcs = fovwcs.deepcopy()

                slicewcs.wcs.ctype = ["LINEAR", "LINEAR",
                                      slicewcs.wcs.ctype[2]]
                slicewcs.wcs.crpix[1] = (ny_slice + 1) / 2
                slicewcs.wcs.crval[1] = (ymin + ymax) / 2 / 3600
                slicewcs.wcs.cdelt[1] = (ymax - ymin) / ny_slice / 3600
                slicewcs_spat = slicewcs.sub(2)

                # World coordinates for the slice
                xworld, yworld = slicewcs_spat.all_pix2world(xslice,
                                                             yslice, 0)
                # FOV pixel coordinates for the slice
                xfov, yfov = fovwcs_spat.all_world2pix(xworld, yworld, 0)

                slicecube = np.zeros((n_z, ny_slice, n_x))
                for islice in range(n_z):
                    ifov = RectBivariateSpline(np.arange(n_y),
                                               np.arange(n_x),
                                               fovcube[islice], kx=1, ky=1)
                    slicecube[islice] = ifov(yfov, xfov, grid=False)

                slicefov = FieldOfView(obj.header,
                                       [obj.meta["wave_min"],
                                        obj.meta["wave_max"]])
                slicefov.detector_header = obj.detector_header
                slicefov.meta["xi_min"] = obj.meta["xi_min"]
                slicefov.meta["xi_max"] = obj.meta["xi_max"]
                slicefov.meta["trace_id"] = sptid
                slicefov.cube = fits.ImageHDU(header=slicewcs.to_header(),
                                              data=slicecube)
                # slicefov.cube.writeto(f"slicefov_{sptid}.fits", overwrite=True)
                slicefov.hdu = spt.map_spectra_to_focal_plane(slicefov)
                if slicefov.hdu is not None:
                    sxmin = slicefov.hdu.header["XMIN"]
                    sxmax = slicefov.hdu.header["XMAX"]
                    symin = slicefov.hdu.header["YMIN"]
                    symax = slicefov.hdu.header["YMAX"]
                    fovimage[symin:symax, sxmin:sxmax] += slicefov.hdu.data

            obj.hdu = fits.ImageHDU(data=fovimage, header=obj.detector_header)

        return obj

    def make_spectral_traces(self):
        """Compute the transformations by interpolation."""
        # nslice = len(self._file["Aperture List"].data)
        # determine echelle order and angle from specified wavelength
        tempres = self._angle_from_lambda()
        self.meta["order"] = tempres["Ord"]
        self.meta["angle"] = tempres["Angle"]

        spec_traces = {}
        for sli in np.arange(self.meta["nslice"]):
            slicename = "Slice " + str(sli + 1)
            spec_traces[slicename] = MetisLMSSpectralTrace(
                self._file,
                spslice=sli, params=self.meta)

        self.spectral_traces = spec_traces

    def rectify_cube(self, hdulist, xi_min=None, xi_max=None, interps=None,
                     **kwargs):
        """
        Rectify an IFU observation into a data cube

        The HDU list (or fits file) must have been created with the
        present OpticalTrain (or an identically configured one).

        Parameters
        ----------
        hdulist : str or fits.HDUList
           an ifu observation created with the present OpticalTrain
        xi_min, xi_max : float [arcsec]
           Spatial limits of the image slicer on the sky. For METIS LMS,
           these values need not be provided by the user.
        interps : list of interpolation functions
           If provided, there must be one for each image extension in `hdulist`.
           The functions go from pixels to the images and can be created with,
           e.g., RectBivariateSpline.
        """
        try:
            inhdul = fits.open(hdulist)
        except TypeError:
            inhdul = hdulist

        # Create interpolation functions
        if interps is None:
            logger.info("Computing interpolation functions")
            interps = make_image_interpolations(inhdul, kx=1, ky=1)

        # Create a common wcs for the rectification
        dwave = from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)
        xi_min = np.min(self.slicelist["left"])
        xi_max = np.max(self.slicelist["right"])
        wave_min = self.meta["wave_min"]
        wave_max = self.meta["wave_max"]
        pixscale = self.meta["pixel_scale"]
        naxis1 = int((xi_max - xi_min) / pixscale) + 1
        naxis2 = len(self.spectral_traces)
        naxis3 = int((wave_max - wave_min)/dwave) + 1
        logger.debug("Cube: %d, %d, %d", naxis1, naxis2, naxis3)
        logger.debug("Xi: %.2f, %.2f", xi_min, xi_max)
        logger.debug("Wavelength: %.3f, %.3f", wave_min, wave_max)
        slicewidth = (self.meta["y_max"] - self.meta["y_min"]) / naxis2

        rectwcs = WCS(naxis=2)
        rectwcs.wcs.ctype = ["WAVE", "LINEAR"]
        rectwcs.wcs.crpix = [1, 1]
        rectwcs.wcs.crval = [wave_min, xi_min]
        rectwcs.wcs.cdelt = [dwave, pixscale]
        rectwcs.wcs.cunit = ["um", "arcsec"]

        cube = np.zeros((naxis3, naxis2, naxis1), dtype=np.float32)
        for i, spt in enumerate(self.spectral_traces.values()):
            spt.wave_min = wave_min
            spt.wave_max = wave_max
            result = spt.rectify(hdulist, interps=interps,
                                 wave_min=wave_min, wave_max=wave_max,
                                 xi_min=xi_min, xi_max=xi_max,
                                 bin_width=dwave)
            cube[:, i, :] = result.data.T

        # FIXME: use wcs object here
        cubehdr = fits.Header()
        cubehdr["INSMODE"] = from_currsys(self.meta["element_name"], self.cmds)
        cubehdr["WAVELEN"] = from_currsys(self.meta["wavelen"], self.cmds)
        cubehdr["CTYPE1"] = "LINEAR"
        cubehdr["CTYPE2"] = "LINEAR"
        cubehdr["CTYPE3"] = "WAVE"
        cubehdr["CRPIX1"] = (naxis1 + 1)/2
        cubehdr["CRPIX2"] = (naxis2 + 1)/2
        cubehdr["CRPIX3"] = 1.
        cubehdr["CRVAL1"] = 0.
        cubehdr["CRVAL2"] = 0.
        cubehdr["CRVAL3"] = self.meta["wave_min"]
        cubehdr["CDELT1"] = pixscale
        cubehdr["CDELT2"] = slicewidth
        cubehdr["CDELT3"] = dwave
        cubehdr["CUNIT1"] = "arcsec"
        cubehdr["CUNIT2"] = "arcsec"
        cubehdr["CUNIT3"] = "um"

        cubehdu = fits.ImageHDU(data=cube, header=cubehdr)
        return cubehdu

    def _angle_from_lambda(self):
        """Determine optimal echelle rotation angle for wavelength."""
        lam = from_currsys(self.meta["wavelen"], self.cmds)
        grat_spacing = self.meta["grat_spacing"]
        wcal = self._file["WCAL"].data
        return echelle_setting(lam, grat_spacing, wcal)


class MetisLMSSpectralTrace(SpectralTrace):
    """SpectralTrace for the METIS LM spectrograph."""

    _class_params = {
        "naxis1": 122,
        "nslice": 28,
        "slicewidth": 0.0207,  # arcsec
        "pixscale": 0.0082,    # arcsec
        "grat_spacing": 18.2,
        "fp2_platescale": 0.303,
    }

    def __init__(self, hdulist, spslice, params, **kwargs):
        polyhdu = hdulist["Polynomial coefficients"]
        params.update(kwargs)
        params["aperture_id"] = spslice
        params["slice"] = spslice
        super().__init__(polyhdu, **params)

        self._file = hdulist
        self.meta["description"] = "Slice " + str(spslice + 1)
        self.meta["trace_id"] = f"Slice {spslice + 1}"
        self.meta.update(params)
        # Provisional:
        self.meta["fov"] = self.fov_grid()

    def fov_grid(self):
        """
        Provide information on the source space volume required by the effect.

        Returns
        -------
        A dictionary with entries `wave_min` and `wave_max`, `x_min`, `y_min`,
        `x_max`, `y_max`. Spatial limits refer to the sky and are given in
        arcsec.
        """
        # TODO: Specify in the warning where the functionality should go!
        warnings.warn("The fov_grid method is deprecated and will be removed "
                      "in a future release. The functionality should be moved"
                      " somewhere else.", DeprecationWarning, stacklevel=2)

        aperture = self._file["Aperture list"].data[self.meta["slice"]]
        x_min = aperture["left"]
        x_max = aperture["right"]
        y_min = aperture["bottom"]
        y_max = aperture["top"]

        filename_det_layout = from_currsys("!DET.layout.file_name", cmds=self.cmds)
        layout = ioascii.read(find_file(filename_det_layout))
        det_lims = {}
        xhw = layout["pixel_size"] * layout["x_size"] / 2
        yhw = layout["pixel_size"] * layout["y_size"] / 2
        det_lims["xd_min"] = min(layout["x_cen"] - xhw)
        det_lims["xd_max"] = max(layout["x_cen"] + xhw)
        det_lims["yd_min"] = min(layout["y_cen"] - yhw)
        det_lims["yd_max"] = max(layout["y_cen"] + yhw)
        wave_min, wave_max = self.get_waverange(det_lims)

        # ..todo: just a hack - xi and x are the same except xi is a quantity
        xi_min = quantify(x_min, u.arcsec)
        xi_max = quantify(x_max, u.arcsec)

        return {"x_min": x_min, "x_max": x_max,
                "y_min": y_min, "y_max": y_max,
                "xi_min": xi_min, "xi_max": xi_max,
                "wave_min": wave_min, "wave_max": wave_max,
                "trace_id": self.trace_id}

    def get_waverange(self, det_mm_lims):
        """Determine wavelength range covered by spec. trace on image plane."""
        xmin = det_mm_lims["xd_min"]
        xmax = det_mm_lims["xd_max"]

        lam0 = from_currsys(self.meta["wavelen"], self.cmds)
        xi0 = 0.
        ymid = self.xilam2y(xi0, lam0)[0]   # estimate y level of trace
        waverange = self.xy2lam(np.array([xmin, xmax]), np.array([ymid, ymid]),
                                grid=False)

        return waverange.min(), waverange.max()

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
        self.xilam2x = Transform2D(matrices["A"].T,
                                   pretransform_x=self.sky2fp,
                                   pretransform_y=self.lam2phase)
        self.xilam2y = Transform2D(matrices["B"].T,
                                   pretransform_x=self.sky2fp,
                                   pretransform_y=self.lam2phase)
        self.xy2lam = Transform2D(matrices["AI"],
                                  posttransform=self.phase2lam)
        self.xy2xi = Transform2D(matrices["BI"],
                                 posttransform=self.fp2sky)

    def get_matrices(self):
        """
        Extract matrix from lms_dist_poly.txt.

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
        spslice = self.meta["slice"]
        order = self.meta["order"]
        angle = self.meta["angle"]
        matnames = ["A", "B", "AI", "BI"]
        matrices = {}

        poly = self.table
        for matid in range(4):
            select = ((poly["Ord"] == order) *
                      (poly["Sli"] == spslice) *
                      (poly["Mat"] == matid))
            if not np.any(select):
                raise KeyError("Combination of Order, Slice not found")

            subpoly = poly[select]
            thematrix = np.zeros((4, 4))
            for i in range(4):
                for j in range(4):
                    sel_ij = (subpoly["Row"] == i) * (subpoly["Col"] == j)
                    thematrix[i, j] = (subpoly["P3"][sel_ij] * angle**3 +
                                       subpoly["P2"][sel_ij] * angle**2 +
                                       subpoly["P1"][sel_ij] * angle +
                                       subpoly["P0"][sel_ij])

            matrices[matnames[matid]] = thematrix

        return matrices

    def lam2phase(self, lam):
        """
        Convert wavelength to phase.

        Phase is lam * order / (2 * grat_spacing).

        Parameters
        ----------
        lam : ndarray (float)
            wavelength (um)

        Returns
        -------
        Phase : ndarray
        """
        return self.meta["order"] * lam / (2 * self.meta["grat_spacing"])

    def phase2lam(self, phase):
        """
        Convert phase to wavelength.

        Wavelength is phase * 2 * grat_spacing / order

        Parameters
        ----------
        phase : ndarray (float)
            phase (dimensionless)

        Returns
        -------
        wavelength : ndarray (um)
        """
        return 2 * self.meta["grat_spacing"] * phase / self.meta["order"]

    def sky2fp(self, xi):
        """Convert position in arcsec to position in FP2."""
        return xi / self.meta["fp2_platescale"]

    def fp2sky(self, fp_x):
        """Convert position in FP2 to position on sky."""
        return fp_x * self.meta["fp2_platescale"]

    def __repr__(self):
        msg = (f"{self.__class__.__name__}({self._file!r}, "
               f"{self.meta['slice']!r}, {self.meta!r})")
        return msg

    def __str__(self):
        msg = (f"<MetisLMSSpectralTrace> \"{self.meta['description']}\" : "
               f"{from_currsys(self.meta['wavelen'], self.cmds)} um : "
               f"Order {self.meta['order']} : Angle {self.meta['angle']}")
        return msg


def echelle_setting(wavelength, grat_spacing, wcal_def):
    """
    Determine optimal echelle rotation angle for wavelength.

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
            wcal = fits.getdata(wcal_def, extname="WCAL")
        except OSError:
            wcal = ioascii.read(wcal_def, comment="^#", format="csv")
    else:
        raise TypeError("wcal_def not in recognised format:", wcal_def)

    # Compute angles, determine which order gives angle closest to zero
    angles = wcal["c0"] * wavelength + wcal["c1"]
    imin = np.argmin(np.abs(angles))

    # Extract parameters
    order = wcal["Ord"][imin]
    angle = angles[imin]

    # Compute the phase corresponding to the wavelength
    phase = wavelength * order / (2 * grat_spacing)

    return {"Ord": order, "Angle": angle, "Phase": phase}


class MetisLMSImageSlicer(ApertureMask):
    """
    Treats the METIS LMS image slicer as an aperture mask effect.

    This helps in building a ``FieldOfView`` object that combines the spatial
    field of the slicer with the spectral range covered by the LMS setting.

    The effect differs from its parent class ``ApertureMask`` in the
    initialisation from the `Aperture List` extension of the trace file
    `!OBS.trace_file`.
    """

    def __init__(self, filename, ext_id="Aperture List", **kwargs):
        filename = find_file(from_currsys(filename, kwargs.get("cmds")))
        ap_hdr = fits.getheader(filename, extname=ext_id)
        ap_list = fits.getdata(filename, extname=ext_id)
        xmin, xmax = ap_list["left"].min(), ap_list["right"].max()
        ymin, ymax = ap_list["bottom"].min(), ap_list["top"].max()
        slicer_dict = {"x": [xmin, xmax, xmax, xmin],
                       "y": [ymin, ymin, ymax, ymax]}
        try:
            kwargs["x_unit"] = ap_hdr["X_UNIT"]
            kwargs["y_unit"] = ap_hdr["Y_UNIT"]
        except KeyError:
            pass

        super().__init__(array_dict=slicer_dict, id="LMS slicer",
                         conserve_image=True, **kwargs)


class MetisLMSEfficiency(TERCurve):
    """
    Computes the grating efficiency (blaze function) for the METIS LMS.

    The procedure is described in E-REP-ATC-MET-1016_1.0. For a given order
    (determined by the central wavelength) the grating efficiency is modelled
    as a squared sinc function of wavelength via the grating angle.
    """

    _class_params = {
        "grat_spacing": 18.2,
        "eff_wid": 0.23,
        "eff_max": 0.75
    }

    def __init__(self, **kwargs):
        # TODO: Refactor these _class_params?
        self.meta = copy.copy(self._class_params)
        assert "grat_spacing" in self.meta, "grat_spacing is missing from self.meta 1"
        super().__init__(**(kwargs | self.meta))
        assert "grat_spacing" in self.meta, "grat_spacing is missing from self.meta 2"

        filename = find_file(self.meta["filename"])
        wcal = fits.getdata(filename, extname="WCAL")
        if "wavelen" in kwargs:
            wavelen = from_currsys(kwargs["wavelen"], kwargs.get("cmds"))
            grat_spacing = self.meta["grat_spacing"]
            ech = echelle_setting(wavelen, grat_spacing, wcal)
            self.meta["order"] = ech["Ord"]
        else:
            wavelen = None

        lam, efficiency = self.make_ter_curve(wcal, wavelen)

        super().__init__(wavelength=lam,
                         transmission=efficiency,
                         emissivity=np.zeros_like(lam),
                         **self.meta)

    def make_ter_curve(self, wcal, wavelen=None):
        """Compute the blaze function for the selected order."""
        order = self.meta["order"]
        eff_wid = self.meta["eff_wid"]
        eff_max = self.meta["eff_max"]

        wcal_ord = wcal[wcal["Ord"] == self.meta["order"]]

        if wavelen is not None:
            lam = np.linspace(wavelen - 0.2, wavelen + 0.2, 1001)
            angle = wcal_ord["c0"] * lam + wcal_ord["c1"]
        else:
            angle = np.linspace(7, -7, 10001)
            lam = wcal_ord["ic0"] * angle + wcal_ord["ic1"]

        phase = order * np.pi * np.sin(np.deg2rad(angle)) * eff_wid
        efficiency = eff_max * np.sinc(phase / np.pi)**2

        return lam, efficiency
