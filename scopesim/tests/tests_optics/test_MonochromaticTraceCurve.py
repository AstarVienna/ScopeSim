import pytest
import numpy as np
from astropy.wcs import WCS
from matplotlib import pyplot as plt

from scopesim.base_classes import PoorMansHeader
from scopesim.optics.monochromatic_trace_curve import MonochromeTraceCurve


PLOTS = False


def _basic_mtc(**kwargs):
    x = np.array([1, 2, 3, 7])
    y = np.array([11, 12, 13, 17])
    s = np.array([-0.5, 0.5, 1.5, 10.])

    wave_min, wave_max = 1.4, 1.5

    return MonochromeTraceCurve(x, y, s, wave_min, wave_max, **kwargs)


@pytest.fixture(scope="function")
def basic_mtc():
    return _basic_mtc()


class TestInit:
    def test_initialises_with_nothing(self):
        with pytest.raises(TypeError):
            MonochromeTraceCurve()

    @pytest.mark.usefixtures("basic_mtc")
    def test_initialises_with_all_coords(self, basic_mtc):
        assert isinstance(basic_mtc, MonochromeTraceCurve)

    def test_initialises_with_non_standard_planes(self):
        mtc = _basic_mtc(aperture_id=2, image_plane_id=1)
        assert mtc.meta["aperture_id"] == 2
        assert mtc.meta["image_plane_id"] == 1


class TestGetHeader:
    @pytest.mark.usefixtures("basic_mtc")
    def test_returns_astropy_header_object_with_basic_input(self, basic_mtc):
        basic_mtc.meta["pixel_size"] = 1    # mm/pix
        print(dict(basic_mtc.header))
        # assert isinstance(basic_mtc.header, fits.Header)
        assert isinstance(basic_mtc.header, PoorMansHeader)

    @pytest.mark.usefixtures("basic_mtc")
    def test_header_reprojects_properly(self):

        if PLOTS:
            xs = [[-1, 0, 1], [2, 3, 4], [-3, -4, -5]]
            ys = [[-1, -1, -1], [0, 0, 0], [1, 1, 1]]
            for x, y, c1, c2 in zip(xs, ys, "rgb", "mck"):
                mtc = MonochromeTraceCurve(x, y, [-1, 0, 1], 1.2, 1.4)

                pixel_size = 0.1
                hdr = mtc.get_header(pixel_size)
                plt.plot(mtc.x, mtc.y, c1)
                plt.plot(mtc.x[0], mtc.y[0], c1+"o")

                xp = [0, hdr["NAXIS1"]]
                yp = [0, 0]
                wcs = WCS(hdr, key="D")
                xw, yw = wcs.all_pix2world(xp, yp, 1)
                plt.plot(xw, yw, c2)
                plt.plot(hdr["CRVAL1D"], hdr["CRVAL2D"], c2+"o")

            plt.show()
