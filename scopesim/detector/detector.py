import numpy as np

from ..base_classes import ImagePlaneBase, DetectorBase
from ..optics import image_plane_utils as imp_utils

from astropy.io import fits


class Detector(DetectorBase):
    def __init__(self, header, **kwargs):
        image = np.zeros((header["NAXIS1"], header["NAXIS2"]))
        self.image_hdu = fits.ImageHDU(header=header, data=image)
        self.meta = {}
        self.meta.update(kwargs)

    def extract(self, image_plane, order=1):
        if not isinstance(image_plane, ImagePlaneBase):
            raise ValueError("image_plane must be an ImagePlane object: {}"
                             "".format(type(image_plane)))

        self.image_hdu = imp_utils.add_imagehdu_to_imagehdu(image_plane.hdu,
                                                            self.hdu, order,
                                                            wcs_suffix="D")

    @property
    def hdu(self):
        self.image_hdu.header.update(self.meta)
        return self.image_hdu

    @property
    def header(self):
        return self.image_hdu.header

    @property
    def data(self):
        return self.image_hdu.data

    @property
    def image(self):
        return self.data
