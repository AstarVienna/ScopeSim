import os
import pytest

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from matplotlib import pyplot as plt

from scopesim.effects.spectral_trace_list import SpectralTraceList
from scopesim.optics.fov_manager import FovVolumeList
from scopesim.tests.mocks.py_objects import trace_list_objects as tlo
from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.base_classes import PoorMansHeader
from scopesim import rc

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SPEC/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]

PLOTS = False


@pytest.fixture(scope="function")
def slit_header():
    return ho._short_micado_slit_header()


@pytest.fixture(scope="function")
def long_slit_header():
    return ho._long_micado_slit_header()


@pytest.fixture(scope="function")
def full_trace_list():
    return tlo.make_trace_hdulist()


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(SpectralTraceList(), SpectralTraceList)

    @pytest.mark.usefixtures("full_trace_list")
    def test_initialises_with_a_hdulist(self, full_trace_list):
        spt = SpectralTraceList(hdulist=full_trace_list)
        assert isinstance(spt, SpectralTraceList)
        assert spt.get_data(2, fits.BinTableHDU)

    def test_initialises_with_filename(self):
        spt = SpectralTraceList(filename="TRACE_MICADO.fits",
                                wave_colname="wavelength", s_colname="xi")
        assert isinstance(spt, SpectralTraceList)


@pytest.mark.skip(reason="Ignoring old Spectroscopy integration tests")
class TestGetFOVHeaders:
    @pytest.mark.usefixtures("full_trace_list", "slit_header")
    def test_gets_the_headers(self, full_trace_list, slit_header):
        spt = SpectralTraceList(hdulist=full_trace_list)
        params = {"pixel_scale": 0.015, "plate_scale": 0.26666,
                  "wave_min": 0.7, "wave_max": 2.5}
        hdrs = spt.get_fov_headers(slit_header, **params)

        # assert all([isinstance(hdr, fits.Header) for hdr in hdrs])
        assert all([isinstance(hdr, PoorMansHeader) for hdr in hdrs])
        # ..todo:: add in some better test of correctness

        if PLOTS:
            # pixel coords
            for hdr in hdrs[::50]:
                xp = [0, hdr["NAXIS1"], hdr["NAXIS1"], 0]
                yp = [0, 0, hdr["NAXIS2"], hdr["NAXIS2"]]
                wcs = WCS(hdr, key="D")
                # world coords
                xw, yw = wcs.all_pix2world(xp, yp, 1)
                plt.fill(xw / hdr["CDELT1D"], yw / hdr["CDELT2D"], alpha=0.2)
            plt.show()

    def test_gets_headers_from_real_file(self):
        slit_hdr = ho._long_micado_slit_header()
        # slit_hdr = ho._short_micado_slit_header()
        wave_min = 1.0
        wave_max = 1.3
        spt = SpectralTraceList(filename="TRACE_15arcsec.fits",
                                s_colname="xi",
                                wave_colname="lam",
                                col_number_start=1,
                                invalid_value=0.,
                                spline_order=1)
        params = {"wave_min": wave_min, "wave_max": wave_max,
                  "pixel_scale": 0.004, "plate_scale": 0.266666667}
        hdrs = spt.get_fov_headers(slit_hdr, **params)
        assert isinstance(spt, SpectralTraceList)

        print(len(hdrs))

        if PLOTS:
            spt.plot(wave_min, wave_max)

            # pixel coords
            for hdr in hdrs[::300]:
                xp = [0, hdr["NAXIS1"], hdr["NAXIS1"], 0]
                yp = [0, 0, hdr["NAXIS2"], hdr["NAXIS2"]]
                wcs = WCS(hdr, key="D")
                # world coords
                xw, yw = wcs.all_pix2world(xp, yp, 1)
                plt.plot(xw, yw, alpha=0.2)
            plt.show()


# class TestApplyTo:
#     def test_fov_setup_base_returns_only_extracted_fov_limits(self):
#         fname = r"F:\Work\irdb\MICADO\TRACE_MICADO.fits"
#         spt = SpectralTraceList(filename=fname, s_colname='xi')
#
#         fvl = FovVolumeList()
#         fvl = spt.apply_to(fvl)
#
#         assert len(fvl) == 17


################################################################################


def test_set_pc_matrix(rotation_ang=0, shear_ang=10):
    n = 100
    im = np.arange(n**2).reshape(n, n)
    hdu = fits.ImageHDU(im)
    hdr_dict = {"CTYPE1": "LINEAR",
                "CTYPE2": "LINEAR",
                "CUNIT1": "deg",
                "CUNIT2": "deg",
                "CDELT1": 1,
                "CDELT2": 1,
                "CRVAL1": 0,
                "CRVAL2": 0,
                "CRPIX1": 0,
                "CRPIX2": 0}
    hdu.header.update(hdr_dict)

    c = np.cos(rotation_ang / 57.29578) * 2
    s = np.sin(rotation_ang / 57.29578) * 2
    t = np.tan(shear_ang / 57.29578)

    n = 5
    pc_dict = {"PC1_1": c + t*s,
               "PC1_2": -s + t*c,
               "PC2_1": s,
               "PC2_2": c}
    det = np.sqrt(np.abs(pc_dict["PC1_1"] * pc_dict["PC2_2"] - \
                         pc_dict["PC1_2"] * pc_dict["PC2_1"]))
    for key in pc_dict:
        pc_dict[key] /= det
    hdu.header.update(pc_dict)
    w = WCS(hdu)

    xd = np.array([0, 10, 10, 0])
    yd = np.array([0, 0, 10, 10])
    xs, ys = w.all_pix2world(xd, yd, 1)

    if PLOTS:
        plt.figure(figsize=(6, 6))
        plt.plot(xd, yd, "o-")
        plt.plot(xs, ys, "o-")
        plt.show()
