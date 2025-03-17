import numpy as np

from astropy.io.fits import ImageHDU
from astropy.wcs import WCS

from ..base_classes import ImagePlaneBase, DetectorBase
from ..optics.image_plane_utils import (add_imagehdu_to_imagehdu,
                                        sky_wcs_from_det_wcs)
from ..utils import get_logger, from_currsys, stringify_dict


logger = get_logger(__name__)


class Detector(DetectorBase):
    def __init__(self, header, cmds=None, **kwargs):
        image = np.zeros((header["NAXIS2"], header["NAXIS1"]))
        self._hdu = ImageHDU(header=header, data=image)
        self.meta = {}
        self.meta.update(header)
        self.meta.update(kwargs)
        self.cmds = cmds

    def extract_from(self, image_plane, spline_order=1, reset=True):
        """Extract HDU from ImagePlane object and add to internal HDU."""
        if reset:
            self.reset()
        if not isinstance(image_plane, ImagePlaneBase):
            raise ValueError("image_plane must be an ImagePlane object, "
                             f"but is: {type(image_plane)}")

        self._hdu = add_imagehdu_to_imagehdu(
            image_plane.hdu, self.hdu, spline_order, wcs_suffix="D")

    def reset(self):
        """Reset internal HDU data to all-zeros-array."""
        # The detector might have been converted to integers by the
        # Quantization effect, so it is not possible to use zeros_like.
        self._hdu.data = np.zeros(self._hdu.data.shape)

    @property
    def hdu(self):
        """Return internal HDU."""
        self._hdu.header.update(stringify_dict(self.meta))

        pixel_scale = from_currsys("!INST.pixel_scale", self.cmds)
        plate_scale = from_currsys("!INST.plate_scale", self.cmds)
        if pixel_scale == 0 or plate_scale == 0:
            logger.warning("Could not create sky WCS.")
        else:
            sky_wcs, _ = sky_wcs_from_det_wcs(
                WCS(self._hdu.header, key="D"), pixel_scale, plate_scale)
            self._hdu.header.update(sky_wcs.to_header())

        return self._hdu

    @property
    def header(self):
        """Return header from internal HDU."""
        return self._hdu.header

    @property
    def data(self):
        """Return data from internal HDU."""
        return self._hdu.data

    @property
    def image(self):
        """Return data from internal HDU."""
        return self.data
