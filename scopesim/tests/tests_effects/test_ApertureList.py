
import numpy as np
from astropy.io import fits

from scopesim.effects import ApertureList, ApertureMask

import matplotlib.pyplot as plt


PLOTS = False


class TestInit:
    def test_initialises_with_nothing(self):
        isinstance(ApertureList(), ApertureList)

    def test_initialises_with_array_dict(self):
        params = {"array_dict": {"id": [0], "left": [-1], "right": [1],
                                 "top": [1], "bottom": [1], "angle": [0],
                                 "conserve_image": [True], "shape": [7]},
                  "x_unit": "arcsec",
                  "y_unit": "arcsec",
                  "pixel_scale": 0.01}
        apl = ApertureList(**params)
        assert isinstance(apl, ApertureList)

    def test_initialises_with_filename(self, mock_path):
        apl = ApertureList(filename=str(mock_path / "test_aperture_list.dat"))
        assert isinstance(apl, ApertureList)


class TestApertures:
    def test_returns_list_of_aperture_masks(self, mock_path):
        apl = ApertureList(filename=str(mock_path / "test_aperture_list.dat"),
                           no_mask=False, pixel_scale=0.01)
        apertures = apl.apertures
        assert all([isinstance(am, ApertureMask) for am in apertures])

        if PLOTS:
            for ii in range(len(apertures)):
                plt.subplot(2, 2, ii+1)
                plt.imshow(apertures[ii].mask.T)
            plt.show()


class TestFovGrid:
    def test_returns_headers(self, mock_path):
        apl = ApertureList(filename=str(mock_path / "test_aperture_list.dat"),
                           no_mask=False, pixel_scale=0.01)
        hdrs = apl.fov_grid()
        assert all([isinstance(hdr, fits.Header) for hdr in hdrs])

        if PLOTS:
            from scopesim.optics.image_plane import ImagePlane
            from scopesim.optics.image_plane_utils import header_from_list_of_xy
            from astropy.io.fits import ImageHDU

            x = np.array([-3, 3]) / 3600.
            y = np.array([-3, 3]) / 3600.
            pixel_scale = 0.01 / 3600
            hdr = header_from_list_of_xy(x, y, pixel_scale)
            implane = ImagePlane(hdr)

            for i in range(4):
                hdu = ImageHDU(data=apl[i].mask.astype(int)+1,
                               header=apl[i].header)
                implane.add(hdu)

            plt.imshow(implane.data, origin="lower")
            plt.show()
