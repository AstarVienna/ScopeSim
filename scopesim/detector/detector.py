from ..optics.image_plane import ImagePlane

from astropy.io import fits


class Detector:
    def __init__(self, header, **kwargs):
        self.image_hdu = fits.ImageHDU(header=header)
        self.meta = kwargs

    def extract(self, image_plane):
        if not isinstance(image_plane, ImagePlane):
            raise ValueError("image_plane must be an ImagePlane object: {}"
                             "".format(type(image_plane)))


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
