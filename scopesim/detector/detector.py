import numpy as np

from ..base_classes import ImagePlaneBase, DetectorBase
from ..optics import image_plane_utils as imp_utils
from ..utils import get_logger, from_currsys, stringify_dict

from astropy.io import fits
from astropy.wcs import WCS


logger = get_logger(__name__)


class Detector(DetectorBase):
    def __init__(self, header, cmds=None, **kwargs):
        image = np.zeros((header["NAXIS2"], header["NAXIS1"]))
        self._hdu = fits.ImageHDU(header=header, data=image)
        self.meta = {}
        self.meta.update(header)
        self.meta.update(kwargs)
        self.cmds = cmds

    def extract_from(self, image_plane, spline_order=1, reset=True):
        if reset:
            self.reset()
        if not isinstance(image_plane, ImagePlaneBase):
            raise ValueError("image_plane must be an ImagePlane object, but is: "
                             f"{type(image_plane)}")

        self._hdu = imp_utils.add_imagehdu_to_imagehdu(image_plane.hdu,
                                                       self.hdu, spline_order,
                                                       wcs_suffix="D")

    def reset(self):
        self._hdu.data = np.zeros_like(self._hdu.data)

    @property
    def hdu(self):
        new_meta = stringify_dict(self.meta)
        self._hdu.header.update(new_meta)

        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        plate_scale = from_currsys("!INST.plate_scale", self.cmds)
        if pixel_scale == 0 or plate_scale == 0:
            logger.warning("Could not create sky WCS.")
        else:
            sky_wcs, _ = imp_utils.sky_wcs_from_det_wcs(
                WCS(self._hdu.header, key="D"), pixel_scale, plate_scale)
            self._hdu.header.update(sky_wcs.to_header())

        return self._hdu

    @property
    def header(self):
        return self._hdu.header

    @property
    def data(self):
        return self._hdu.data

    @property
    def image(self):
        return self.data
