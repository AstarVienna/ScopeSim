import numpy as np
import pytest
from pytest import approx

from astropy import units as u

from scopesim.optics.fov import FieldOfView
from scopesim.optics.image_plane import ImagePlane

from scopesim.tests.mocks.py_objects import source_objects as src
from scopesim.tests.mocks.py_objects import header_objects as hdrs

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


@pytest.fixture(scope="function")
def comb_src():
    return src._combined_source()


@pytest.fixture(scope="function")
def fov_hdr():
    return hdrs._fov_header()


@pytest.fixture(scope="function")
def implane_hdr():
    return hdrs._implane_header()


@pytest.mark.usefixtures("comb_src", "fov_hdr", "implane_hdr")
class TestInteractionBetweenSourceFOVImagePlane:
    def test_can_extract_the_source_in_a_fov(self, fov_hdr, comb_src,
                                             implane_hdr):

        from scipy.ndimage.interpolation import zoom
        comb_src.fields[3].data = zoom(comb_src.fields[3].data, (1.5, 1), order=1)
        print(comb_src.fields[3].data.shape)
        print(dict(comb_src.fields[3].header))

        # angle = 30
        # deg2rad = 3.1415926 / 180
        # comb_src.fields[3].header["PC1_1"] = np.cos(angle * deg2rad)
        # comb_src.fields[3].header["PC1_2"] = np.sin(angle * deg2rad)
        # comb_src.fields[3].header["PC2_1"] = -np.sin(angle * deg2rad)
        # comb_src.fields[3].header["PC2_2"] = np.cos(angle * deg2rad)

        as2deg = u.arcsec.to(u.deg)

        fov = FieldOfView(fov_hdr, waverange=[0.5, 2.5]*u.um)
        fov.hdu.header["CRVAL1"] += 1.5 * as2deg
        fov.hdu.header["CRVAL2"] -= 1.5 * as2deg
        fov.hdu.header["CRVAL1D"] += 30
        fov.hdu.header["CRVAL2D"] -= 30

        imp = ImagePlane(implane_hdr)
        fov.extract_from(comb_src)
        fov.view()
        imp.add(fov.hdu, wcs_suffix="D")

        ipt = np.sum(fov.fields[0]["flux"]) + np.sum(fov.fields[1].data)
        opt = np.sum(imp.image)

        assert ipt == approx(opt)

        if PLOTS:
            plt.subplot(131)
            plt.imshow(comb_src.fields[3].data.T, origin="lower", norm=LogNorm())

            plt.subplot(132)
            plt.imshow(fov.image.T, origin="lower", norm=LogNorm())

            plt.subplot(133)
            plt.imshow(imp.image.T, origin="lower", norm=LogNorm())
            plt.show()
