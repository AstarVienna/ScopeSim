from copy import deepcopy
import numpy as np

from ..base_classes import ImagePlaneBase, DetectorBase
from ..optics import image_plane_utils as imp_utils
from .. import utils

from astropy.io import fits


class Detector(DetectorBase):
    def __init__(self, header, **kwargs):
        image = np.zeros((header["NAXIS2"], header["NAXIS1"]))
        self._hdu = fits.ImageHDU(header=header, data=image)
        self.meta = {}
        self.meta.update(header)
        self.meta.update(kwargs)

    def extract_from(self, image_plane, spline_order=1, reset=True):
        if reset:
            self.reset()
        if not isinstance(image_plane, ImagePlaneBase):
            raise ValueError("image_plane must be an ImagePlane object: {}"
                             "".format(type(image_plane)))

        self._hdu = imp_utils.add_imagehdu_to_imagehdu(image_plane.hdu,
                                                       self.hdu, spline_order,
                                                       wcs_suffix="D")

    def reset(self):
        self._hdu.data = np.zeros(self._hdu.data.shape)

    @property
    def hdu(self):
        new_meta = utils.stringify_dict(self.meta)
        self._hdu.header.update(new_meta)
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
