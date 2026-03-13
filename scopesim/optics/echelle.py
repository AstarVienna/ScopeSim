import dataclasses

import numpy as np
import astropy.units as u

from ..utils import get_logger
logger = get_logger(__name__)

def spectrograph_factory(min_wave: float|u.Quantity, max_wave: float|u.Quantity, focal_len: float|u.Quantity,
                         design_res: float, echelle_angle: float|u.Quantity, min_order: int, max_order: int,
                         echelle_groove_length: float|u.Quantity,
                         pix_per_res_elem: float, disp_npix: int, xdisp_npix: int, pix_size: float|u.Quantity,
                         xdisp_groove_length: float|u.Quantity = 0.0, xdisp_beta_center: float|u.Quantity = 0.0):
    """
    
    Parameters
    ----------
    min_wave: float|u.Quantity
        Minimum wavelength. If float, assumed to be in nm.
    max_wave: float|u.Quantity
        Maximum wavelength. If float, assumed to be in nm.
    focal_len: float|u.Quantity
        Focal length. If float, assumed to be in mm.
    design_res: float|int
        Design resolution (dimensionless).
    echelle_angle: float|u.Quantity
        Echelle angle. If float, assumed to be in degrees.
    min_order: int
        Minimum order.
    max_order: int
        Maximum order.
    echelle_groove_length: float|u.Quantity
        Echelle groove length. If float, assumed to be in mm.
    pix_per_res_elem: float|int
        Pixels per resolution element.
    disp_npix: int
        Number of dispersion pixels.
    xdisp_npix: int
        Number of cross-dispersion pixels.
    pix_size: float|u.Quantity
        Pixel size. If float, assumed to be in um.
    xdisp_groove_length: float|u.Quantity
        Cross disperser groove length. If float, assumed to be in mm.
    xdisp_beta_center: float|u.Quantity
        Cross disperser beta center angle. If float, assumed to be in degrees.

    Returns
    -------

    """
    # Convert primitive types to quantities with appropriate units
    if not isinstance(min_wave, u.Quantity):
        min_wave = min_wave * u.nm
    if not isinstance(max_wave, u.Quantity):
        max_wave = max_wave * u.nm
    if not isinstance(focal_len, u.Quantity):
        focal_len = focal_len * u.mm
    if not isinstance(echelle_angle, u.Quantity):
        echelle_angle = echelle_angle * u.deg
    if not isinstance(echelle_groove_length, u.Quantity):
        echelle_groove_length = echelle_groove_length * u.mm
    if not isinstance(pix_size, u.Quantity):
        pix_size = pix_size * u.mm
    if not isinstance(xdisp_groove_length, u.Quantity):
        xdisp_groove_length = xdisp_groove_length * u.mm
    if not isinstance(xdisp_beta_center, u.Quantity):
        xdisp_beta_center = xdisp_beta_center * u.deg

    x_disp_len = (xdisp_npix*pix_size).to(u.mm)

    if not np.isfinite(xdisp_groove_length):
        nu, xsdisp_angle = GratingSetup.estimate_xdisp_freq_and_angle(min_wave, max_wave, x_disp_len, focal_len)
        cross_disperser = GratingSetup(alpha=xsdisp_angle, beta_center=xsdisp_angle, groove_length=1 / nu,
                                       grating_type='vph', vph_center_wave=(min_wave + max_wave) / 2)
    else:
        if xdisp_beta_center != 0:  # TODO probably hsould be a nan, but use of 0 angle incidence isn't likely
            cross_disperser = GratingSetup(alpha=xdisp_beta_center, beta_center=xdisp_beta_center,
                                           groove_length=xdisp_groove_length, grating_type='vph',
                                           vph_center_wave=xdisp_groove_length.to(u.nm) * 2 * np.sin(xdisp_beta_center))
        else:
            cross_disperser = GratingSetup(groove_length=xdisp_groove_length,
                                           guess_littrow=(min_wave, max_wave, x_disp_len, focal_len),
                                           grating_type='vph', vph_center_wave=(min_wave + max_wave) / 2)

    echelle_grating = GratingSetup(alpha=echelle_angle, beta_center=echelle_angle, delta=echelle_angle,
                                   groove_length=echelle_groove_length)

    return SpectrographSetup(order_range=(min_order, max_order),
                                   design_res=design_res,
                                   pixels_per_res_elem=pix_per_res_elem,
                                   focal_length=focal_len,
                                   grating=echelle_grating,
                                   detector=Detector(disp_npix, xdisp_npix, pix_size),
                                   cross_disperser=cross_disperser)


class GratingSetup:
    def __init__(
            self, *,
            alpha: u.Quantity = None,
            delta: u.Quantity = None,
            beta_center: u.Quantity = None,
            groove_length: u.Quantity = None,
            empiric_efficiency_factor: float = 1.0,
            guess_littrow: tuple[u.Quantity,u.Quantity,u.Quantity,u.Quantity] = None,
            grating_type: str = 'echelle',
            vph_center_wave: u.Quantity = None
    ):
        """
        Simulation of an echelle/echellete grating for spectrometer.

        :param float alpha: incident angle in radians
        :param float delta: blaze angle in radians
        :param float beta_center: reflectance angle in radians
        :param u.Quantity groove_length: d, length/groove in units of wavelength (same as schroeder sigma)
        :param float empiric_efficiency_factor: empirically determined blaze factor to correct for grating efficiency relative to ideal peak
        :param bool guess_littrow: if not None, passed onto estimate_xdisp_angle_with_groove_len
        """
        self.grating_type = grating_type.lower()
        assert self.grating_type in ('echelle', 'vph'), f'grating_type must be echelle or vph, not {self.grating_type}'

        if guess_littrow:
            alpha, beta_center = self.estimate_xdisp_angle_with_groove_len(guess_littrow[0], guess_littrow[1], guess_littrow[2],
                                              guess_littrow[3], groove_length)
            self.alpha = alpha
            self.delta = None if self.grating_type == 'vph' else beta_center
            self.beta_center = beta_center
        else:
            self.alpha = alpha.to(u.rad)
            self.delta = delta.to(u.rad) if delta is not None else None
            self.beta_center = beta_center.to(u.rad)
        self.d = groove_length.to(u.mm)
        self.empiric_efficiency_factor = empiric_efficiency_factor

        # Approximate reality vph_constant = np.pi * deltan * d / bragg_angle
        # Empircally we would use with peak efficiency near center of range
        self.vph_constant = np.pi/2*vph_center_wave.to(u.nm) if self.grating_type == 'vph' else None

    @staticmethod
    def estimate_xdisp_angle_with_groove_len(
            l_min: float | u.Quantity,
            l_max: float | u.Quantity,
            detector_length: float | u.Quantity,
            focal_length: float | u.Quantity,
            d: u.Quantity) -> tuple[u.Quantity, u.Quantity]:
        """
        Estimate cross-disperser angle and beta center for a given groove length.
        
        This method calculates the optimal incident angle (alpha) and central reflection 
        angle (beta) for a cross-disperser grating to match the wavelength range to the 
        detector length at the focal plane.
        
        Parameters
        ----------
        l_min : float | u.Quantity
            Minimum wavelength. If float, assumed to be in nm.
        l_max : float | u.Quantity
            Maximum wavelength. If float, assumed to be in nm.
        detector_length : float | u.Quantity
            Length of the detector. If float, assumed to be in mm.
        focal_length : float | u.Quantity
            Focal length of the system. If float, assumed to be in mm.
        d : u.Quantity
            Groove length (spacing) of the grating.
            
        Returns
        -------
        tuple[u.Quantity, u.Quantity]
            A tuple containing:
            - alpha: incident angle in radians
            - beta_center: central reflection angle in radians
        """
        # Convert floats to quantities with appropriate units
        if not isinstance(l_min, u.Quantity):
            l_min = l_min * u.nm
        if not isinstance(l_max, u.Quantity):
            l_max = l_max * u.nm
        if not isinstance(detector_length, u.Quantity):
            detector_length = detector_length * u.mm
        if not isinstance(focal_length, u.Quantity):
            focal_length = focal_length * u.mm
        if not isinstance(d, u.Quantity):
            d = d * u.mm

        # assume we want a dispersion at the midpoint that matches the range to the detector length
        t = detector_length / focal_length
        k = (l_max - l_min) / d
        x = k * (1 + np.sqrt(1 + t ** 2)) / 2 / t
        u = np.arccos(np.sqrt(x))
        v = np.arcsin(k / 2 / np.sqrt(x))
        beta_max = u+v
        beta_min = u-v
        alpha = np.arcsin(l_max / d - np.sin(beta_max))
        return alpha.to(u.rad), ((beta_min + beta_max) / 2).to(u.rad)

    @staticmethod
    def estimate_xdisp_freq_and_angle(l_min: float | u.Quantity,
                                     l_max: float | u.Quantity,
                                     detector_length: float | u.Quantity,
                                     focal_length: float | u.Quantity) -> tuple[u.Quantity, u.Quantity]:
        """
        Estimate cross-disperser littrow angle and groove frequency.

        This method calculates the analytical solution to the system of equations

        nu l_max == Sin(beta)+Sin(beta+deltabeta) nu l_min == Sin(beta)+Sin(beta-deltabeta)
        cos[deltabeta] = f/sqrt(f^2 + (d/2)^)
        sin[deltabeta] = (d/2)/sqrt(f^2 + (d/2)^)

        to match the initial and final wavelengths to the extent of the detector length


        Parameters
        ----------
        l_min : float | u.Quantity
            Minimum wavelength. If float, assumed to be in nm.
        l_max : float | u.Quantity
            Maximum wavelength. If float, assumed to be in nm.
        detector_length : float | u.Quantity
            Length of the detector. If float, assumed to be in mm.
        focal_length : float | u.Quantity
            Focal length of the system. If float, assumed to be in mm.

        Returns
        -------
        tuple[u.Quantity, u.Quantity]
            A tuple containing:
            - ruling frequency
            - littrow angle
        """
        if not isinstance(l_min, u.Quantity):
            l_min = l_min * u.nm
        if not isinstance(l_max, u.Quantity):
            l_max = l_max * u.nm
        if not isinstance(detector_length, u.Quantity):
            detector_length = detector_length * u.mm
        if not isinstance(focal_length, u.Quantity):
            focal_length = focal_length * u.mm

        f = focal_length
        p = detector_length
        li, lf = l_min, l_max
        a = 4*f**2 + p**2
        rta = np.sqrt(a)
        b = lf - li
        bsq = b**2
        c = lf**2 + li**2
        d = 4 * f**2 * bsq + 2 * f * rta * bsq + p**2 * c
        beta = np.arccos(((2 * f + rta) * b) / np.sqrt(2 * d))
        nu = (np.sqrt(2) * p * (p**2 + 4 * f * (2 * f + rta))) / ((2 * f + rta) * np.sqrt(a * d))
        return nu.to('mm-1'), beta.to('rad')

    def __str__(self) -> str:
        return (f"alpha={np.rad2deg(self.alpha):.2f}\n"
                f"delta={np.rad2deg(self.delta):.2f}\n"
                f"beta={np.rad2deg(self.beta_center):.2f}\n"
                f"l={self.d.to('mm'):.2f}/l ({1 / self.d.to('mm'):.2f})")

    def blaze(self, beta, m):
        """
        Blaze throughput function follows Casini & Nelson 2014 J Opt Soc Am A eq 25 with notation modified to
        match Schroder.
        :param beta: reflectance angle in radians
        :param m: order
        :return: grating throughput out of 1
        """
        if self.grating_type == 'vph':
            raise ValueError('Cannot calculate blaze for vph gratings.')
        k = np.cos(beta) * np.cos(self.alpha - self.delta) / (np.cos(self.alpha) * np.cos(beta - self.delta))
        k[k > 1] = 1
        # q1 = np.cos(alpha) / np.cos(alpha - delta)
        # q3=np.cos(delta)-np.sin(delta)*np.cot((alpha+beta)/2)
        # return k*np.sinc(m*q1*q3)**2
        q4 = np.cos(self.delta) - np.sin(self.delta) / np.tan((self.alpha + beta) / 2)
        rho = np.cos(self.delta) if self.alpha < self.delta else np.cos(self.alpha) / np.cos(self.alpha - self.delta)
        # 2 different rho depending on whether alpha or delta is larger
        ret = k * np.sinc((m * rho * q4).value) ** 2  # omit np.pi as np.sinc includes it
        return self.empiric_efficiency_factor * ret

    def vph_efficiency(self, wavelength):
        """
        :param wavelength: wavelength(s) as u.Quantity units of length
        :return: efficiency of the grating
        """
        if self.grating_type != 'vph':
            raise ValueError('Cannot calculate vph_effiency for echelle.')

        return self.empiric_efficiency_factor * np.sin((self.vph_constant/wavelength).decompose().value)**2

    def beta(self, wave, m):
        """
        :param wave: wavelength(s) as u.Quantity units of length
        :param m: order
        :return: reflectance angle in radians
        """
        return np.arcsin(m * wave / self.d - np.sin(self.alpha))

    def wave(self, beta, m):
        """
        :param beta: reflectance angle in radians
        :param m: order
        :return: the wavelength of beta in that order
        """
        return (self.d * (np.sin(beta) + np.sin(self.alpha)) / m).to(u.nm)

    def limiting_resolution(self, beam_diameter: u.Quantity, order):
        """
        :param u.Quantity beam_diameter: diameter of the incoming beam
        :param order: order
        :return: the limiting resolution of the grating configuration
        """
        return order * beam_diameter / (self.d * np.cos(self.alpha))

    def effective_resolution(self,
                       beam_diameter: u.Quantity,
                       order,
                       wave: u.Quantity,
                       phi: float,
                       telescope_diameter: u.Quantity
                       ):
        """
        :param u.Quantity beam_diameter: diameter of the incoming beam
        :param order: order
        :param wave: wavelength(s) as u.Quantity
        :param phi: angular slit width (small angle approx: width/tele_f_len) in radians
        :param u.Quantity telescope_diameter: telescope diameter
        :return: the effective resolution of the grating configuration
        """
        return self.limiting_resolution(beam_diameter, order) * wave / (phi * telescope_diameter)

    def angular_dispersion(self, m, beta):
        """
        :param m: order
        :param beta: reflectance angle
        :return: angular dispersion [rad/wavelength], Schroder A dbeta/dlambda
        """
        return m / (self.d * np.cos(beta)) * u.rad



@dataclasses.dataclass
class Detector:
    n_pix_x: int
    n_pix_y: int
    pixel_size: u.Quantity

class SpectrographSetup:
    def __init__(
            self,
            order_range: tuple,
            design_res: float,
            pixels_per_res_elem: float,
            focal_length: u.Quantity,
            grating: GratingSetup,
            detector: Detector,
            cross_disperser: GratingSetup = None,
    ):
        """

        We assume that the system is designed to place pixels_per_res_elem pixels across the width of the dispersed
        slit image at some fiducial wavelength and that the resolution element has some known width there
        (and that the slit image is gaussian)

        :param tuple order_range: order range of the spectrograph
        # :param u.Quantity final_wave: longest wavelength at the edge of detector
        :param float pixels_per_res_elem: number of pixels per resolution element of spectrometer
        :param u.Quantity focal_length: the focal length of the detector
        :param GratingSetup grating: configured grating
        :param MKIDDetector detector: configured detector
        :return spectrograph simulation object
        """
        assert len(order_range) == 2, 'order_range must be a tuple of length 2.'
        if order_range[0] > order_range[1]:
            order_range = order_range[::-1]
        self.m0 = order_range[0]
        self.m_max = order_range[1]
        self.l0 = grating.wave(grating.beta_center, order_range[0])*(1/order_range[0]/2+1)
        self.design_res = design_res
        self.grating = grating
        self.detector = detector
        self.focal_length = focal_length.to(u.mm)
        self.pixel_scale = np.arctan(self.detector.pixel_size / self.focal_length)
        self.beta_central_pixel = self.grating.beta_center
        self.nord = int(self.m_max - self.m0 + 1)
        self.nominal_pixels_per_res_elem = pixels_per_res_elem
        self.cross_disperser = cross_disperser
        self._orders = None

        self.nondimensional_lsf_width = 1 / self.design_res
        logger.debug(f'\nThe spectrograph has been setup with the following properties:'
                     f'\n\tl0: {self.l0}'
                     # f'\n\tR0: {self.detector.design_R0}'
                     f'\n\tOrders: {self.orders}'
                     f'\n\tFocal length: {self.focal_length}'
                     f'\n\tIncidence angle: {np.rad2deg(self.grating.alpha):.3f}'
                     f'\n\tReflectance angle: {np.rad2deg(self.beta_central_pixel):.2f}\n'
                     f'\n\tGroove length: {self.grating.d:.2f}'
                     f'\n\t# of pixels: {self.detector.n_pix_x}x{self.detector.n_pix_y}'
                     f'\n\tPixel size: {self.detector.pixel_size}'
                     f'\n\tPixels per res. element: {self.nominal_pixels_per_res_elem}')

    def set_beta_center(self, beta, littrow: bool = False):
        """
        :param beta: reflectance angle in degrees
        :param littrow: whether alpha=beta
        :return: changes the reflectance angle of the central pixel (may change alpha if littrow)
        """
        if not isinstance(beta, u.Quantity):
            beta *= u.deg
        self.beta_central_pixel = beta
        if littrow:
            self.grating.alpha = beta

    @property
    def orders(self):
        try:
            assert self._orders[0] == (self.m0, self.m_max)
        except (TypeError, AssertionError):
            self._orders = (self.m0, self.m_max), np.arange(self.m0, self.m_max + 1, dtype=int)
        return self._orders[1]

    @property
    def minimum_wave(self):
        return self.central_wave(self.m_max) - self.fsr(self.m_max) / 2

    def info_str(self):
        gstr = str(self.grating)
        betactr = np.rad2deg(self.grating.beta(self.central_wave(self.m0), self.m0))
        gstr += f'\nb={betactr:.2f}'
        ret = [f'    {x}' for x in gstr.split('\n')]
        ret.insert(0, 'Grating:')

        beta_ext = np.rad2deg(self.beta_for_x_pixel(np.array([0, self.detector.n_pix_x-1]) + .5))
        ret.append("    beta extent: {:.2f} - {:.2f}".format(*beta_ext))
        ret.append('Orders:')
        for o in self.orders[::-1]:
            w_c = self.central_wave(o)
            w_i = w_c - self.fsr(o) / 2
            w_f = w_c + self.fsr(o) / 2
            p_i = self.wavelength_to_x_pixel(w_i, o)
            p_f = self.wavelength_to_x_pixel(w_f, o)
            ret.append(f"    m{o:2} @ {w_c:.0f}: {w_i:.0f} - {w_f:.0f}, {p_i:.0f} - {p_f:.0f}")
        return ret

    def x_pixel_for_beta(self, beta):
        """
        :param beta: reflectance angle in radians
        :return: pixel at beta
        """
        delta_angle = np.tan(beta - self.beta_central_pixel)
        return self.focal_length * delta_angle / self.detector.pixel_size + self.detector.n_pix_x / 2

    def beta_for_x_pixel(self, pixel):
        """
        :param pixel: pixel index
        :return: reflectance angle (radians) at pixel
        """
        center_offset = self.detector.pixel_size * (pixel - self.detector.n_pix_x / 2)
        return self.beta_central_pixel + np.arctan(center_offset / self.focal_length)

    def y_pixel_for_beta(self, beta):
        """
        :param beta: reflectance angle in radians
        :return: pixel at beta
        """
        delta_angle = np.tan(beta - self.cross_disperser.beta_center)
        return self.focal_length * delta_angle / self.detector.pixel_size + self.detector.n_pix_y / 2

    def beta_for_y_pixel(self, pixel):
        """
        :param pixel: pixel index
        :return: reflectance angle (radians) at pixel
        """
        center_offset = self.detector.pixel_size * (pixel - self.detector.n_pix_y / 2)
        return self.cross_disperser.beta_center + np.arctan(center_offset / self.focal_length)

    @property
    def max_beta_m0(self):
        """
        :return: largest reflectance angle (which is at the initial order)
        """
        return self.grating.beta(self.l0, self.m0)

    @property
    def min_beta_mmax(self):
        """
        :return: smallest reflectance angle (which is at the final order)
        """
        return self.grating.beta(self.minimum_wave, self.m_max)

    def blaze(self, wave):
        """
        :param wave: wavelength
        :return: blaze throughput out of 1
        """
        return self.grating.blaze(self.grating.beta(wave, self.orders[:, None]), self.orders[:, None])

    def xdisp_efficiency(self, wave):
        """
        :param wave: wavelength
        :return: efficiency of the cross disperser
        """
        return self.cross_disperser.vph_efficiency(wave)

    def mean_blaze_eff_est(self, n=10):
        """
        :param n:
        :return:
        """
        edges = self.edge_wave(fsr=True)
        detector_edges = self.edge_wave(fsr=False)

        ow = np.array([np.select([detector_edges > edges, detector_edges <= edges],
                                 [detector_edges, edges])[:, 0],
                       np.select([detector_edges < edges, detector_edges > edges],
                                 [detector_edges, edges])[:, 1]]).T * u.nm
        v = self.blaze(np.array(list(map(lambda x: np.linspace(*x, num=n), ow))) * u.nm).mean(1)
        return v.value

    def order_mask(self, wave, fsr_edge: bool = False):
        """
        :param wave: wavelength(s) as u.Quantity
        :param bool fsr_edge: True to mask at the FSR, goes to detector edge if not
        :return: a boolean array [norders, wave.size] where true means wavelengths are in that order
        """
        if fsr_edge:
            o = self.orders[:, None]
            c_wave = self.x_pixel_to_wavelength(self.detector.n_pix_x / 2, o)
            fsr = c_wave / o
            return np.abs(wave - c_wave) < fsr / 2
        else:
            x = self.wavelength_to_x_pixel(wave, self.orders[:, None])
            return (x >= 0) & (x < self.detector.n_pix_x)

    def edge_wave(self, fsr=True):
        """
        :param fsr: True to return the FSR edge wavelengths, False to return the detector edges
        :return: wavelengths at the edges of the detector (or FSR) for each order, always include central pixel wavelength
        """
        pix = np.array([0, self.detector.n_pix_x // 2, self.detector.n_pix_x-1]) + .5
        fiducial_waves = self.x_pixel_to_wavelength(pix, self.orders[:, None])
        if not fsr:
            return fiducial_waves[:, [0, -1]]

        central_fsr = fiducial_waves[:, 1] / self.orders
        fsr_edges = (u.Quantity([-central_fsr / 2, central_fsr / 2]) + fiducial_waves[:, 1]).T

        return fsr_edges

    def central_wave(self, order):
        """
        :param order: order number
        :return: wavelength at the center of the order
        """
        l0_center = self.l0 / (1 + 1 / (2 * self.m0))
        return l0_center * self.m0 / order

    def fsr(self, order):
        """
        :param order: order number
        :return: free spectral range of the order
        """
        return self.central_wave(order) / order

    def wavelength_to_x_pixel(self, wave, m):
        """
        :param wave: wavelength(s) as u.Quantity
        :param m: order
        :return: fractional pixel location of given wavelength
        """
        return self.x_pixel_for_beta(self.grating.beta(wave, m))

    def x_pixel_to_wavelength(self, pixel, m):
        """
        :param pixel: pixel
        :param m: order
        :return: wavelength for given pixel as u.Quantity
        """
        return self.grating.wave(self.beta_for_x_pixel(pixel), m)

    def x_pixel_wavelengths(self, pixel_phase=.5):
        """
        :param pixel_phase: phase of the pixel (0-1), default center (.5)
        :return:  array of pixel center wavelengths for every order

        Note that wavelengths will be computed outside each order's FSR.
        NB this is approximately equal to the simple approximation:
        (np.linspace(-.5, .5, num=n_pixels) * self.fsr(m0) + self.central_wave(m0)) * (m0/self.orders)[:,None]
        """
        assert 0<=pixel_phase<=1, 'pixel_phase must be between 0 and 1'
        return self.x_pixel_to_wavelength(np.arange(self.detector.n_pix_x) + pixel_phase, self.orders[:, None])

    def wavelength_to_y_pixel(self, wave):
        return self.y_pixel_for_beta(self.cross_disperser.beta(wave, 1))

    def y_pixel_to_wavelength(self, pixel, m):
        """
        :param pixel: pixel
        :param m: order
        :return: wavelength for given pixel as u.Quantity
        """
        return self.cross_disperser.wave(self.beta_for_y_pixel(pixel), 1)

    @property
    def dl_pix_max_wave(self):
        """
        :return: maximum change in wavelength in any pixel
        """
        return self.pixel_scale / self.grating.angular_dispersion(self.m0, self.max_beta_m0)

    @property
    def dl_pix_min_wave(self):
        """
        :return: minimum change in wavelength in any pixel
        """
        return self.pixel_scale / self.grating.angular_dispersion(self.m_max, self.min_beta_mmax)

    @property
    def dl_pixel(self):
        """
        :return: change in wavelength for every pixel
        """
        return self.pixel_scale / self.angular_dispersion

    # def pixel_rescale(self, oversampling):
    #     """
    #     :param oversampling: factor by which to oversample smallest wavelength extent
    #     :return: The sample size in wavelength units for every pixel. Every MKID resolution width divided by the
    #              total # of samples that are in the largest width, which is the largest width divided by the
    #              smallest sample size.
    #     """
    #     return self.dl_mkid_pixel * self.sampling(oversampling) / self.dl_mkid_max
    #
    # def pixel_samples_frac(self, oversampling):
    #     """
    #     :return: number of samples for every pixel, retrieved by dividing change in wavelength for a pixel
    #              by the sample size for that pixel.
    #     """
    #     return (self.dl_pixel / self.pixel_rescale(oversampling)).si.value
    #
    # def pixel_max_npoints(self, oversampling):
    #     """
    #     :return: the maximum number of points in any given pixel, as an integer value, ensuring there is at least one
    #              sample at pixel center
    #     """
    #     pixel_max_npoints = np.ceil(self.pixel_samples_frac(oversampling).max()).astype(int)
    #     if not pixel_max_npoints % 2:  # ensure there is a point at the pixel center
    #         pixel_max_npoints += 1
    #     return pixel_max_npoints

    @property
    def angular_cross_dispersion(self):
        """
        :return: The angular dispersion at the center of each pixel for cross dispersion
        """
        beta = self.beta_for_y_pixel(np.arange(self.detector.n_pix_y) + .5)
        return self.cross_disperser.angular_dispersion(1, beta)

    @property
    def angular_dispersion(self):
        """
        :return: The angular dispersion at the center of each pixel for each order (nord, npixel)
        """
        beta = self.beta_for_x_pixel(np.arange(self.detector.n_pix_x) + .5)
        return self.grating.angular_dispersion(self.orders[:, None], beta)

    @property
    def implied_design_res(self):
        """
        :return: design resolution for spectrometer
        Assumes that m0 (the longest wavelength order) FSR fills detector with some sampling
        """
        dlambda = self.fsr(self.m0) / self.detector.n_pix_x * self.nominal_pixels_per_res_elem
        return self.l0 / dlambda

    @property
    def average_res(self):
        """
        :return:
        """
        w = self.edge_wave(fsr=False)
        return (w.mean(1) / (np.diff(w, axis=1).T / self.detector.n_pixels * self.nominal_pixels_per_res_elem)).ravel()

    def plot_echellogram(self, center_orders=True, title='', blaze=False, cross_disperse=False):
        if cross_disperse and self.cross_disperser is None:
            raise ValueError('No cross disperser defined')
        if cross_disperse:
            raise NotImplementedError('Cross dispersion not yet implemented')

        import matplotlib.pyplot as plt
        w = self.x_pixel_wavelengths()
        b = self.blaze(w)

        fig, axes = plt.subplots(2 + int(blaze), 1, figsize=(6, 10 + 4 * int(blaze)))
        if title:
            plt.suptitle(title)
        plt.sca(axes[0])
        plt.title(f'a={self.beta_central_pixel:.1f} m={self.m0}-{self.m_max}')
        fsr_edges = self.edge_wave(fsr=True)
        for ii, i in enumerate(self.orders):
            waves = w[ii, [0, self.detector.n_pix_x // 2, -1]]
            plt.plot(self.wavelength_to_x_pixel(waves, i), [i] * 3, '*', color=f'C{ii}')
            plt.plot(self.wavelength_to_x_pixel(fsr_edges[ii], i), [i] * 2, '.', color=f'C{ii}')
        plt.xlabel('Pixel')
        plt.ylabel('Order')
        plt.sca(axes[1])
        for ii, i in enumerate(self.orders):
            waves = w[ii, [0, self.detector.n_pix_x // 2, -1]]
            oset = waves[1] if center_orders else 0
            plt.plot(waves - oset, [i] * 3, '*', color=f'C{ii}')
            plt.plot(fsr_edges[ii] - oset, [i] * 2, '.', color=f'C{ii}',
                     label=f'$\\lambda=${waves[1]:.0f}')
        plt.legend()
        plt.xlabel('Center relative wavelength (nm)')
        plt.ylabel('Order')
        if blaze:
            plt.sca(axes[2])
            plt.plot(w.T, b.T)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Blaze')
        plt.tight_layout()
        plt.show()


# GRATING_CATALOG = np.loadtxt('newport_masters.txt', delimiter=',',
#                              dtype=[('name', 'U10'), ('l', 'f4'), ('blaze', 'f4'),
#                                     ('width', 'f4'), ('height', 'f4'), ('stock', 'U10')])
# GRATING_CATALOG['l'] = 1e6/GRATING_CATALOG['l']
#
# NEWPORT_GRATINGS = {x['name']: GratingSetup(
#     0,
#     (x['blaze']*u.deg).to(u.rad).value,
#     0,
#     x['l'] * u.nm) for x in GRATING_CATALOG}


# ss = SpectrographSetup((18,29),
#                   515*u.nm,
#                   4.7,
#                   225*u.mm,
#                   GratingSetup(alpha=np.deg2rad(64.2), beta_center=np.deg2rad(64.2), delta=np.deg2rad(64.2), groove_length=u.mm/200),
#                   Detector(4096,4096,15*u.micron),
#                   cross_disperser=GratingSetup(groove_length=u.mm/1000, guess_littrow=(310*u.nm, 515*u.nm, (4096-20)*0.015*u.mm, 225*u.mm)),
#                   )