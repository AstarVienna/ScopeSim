"""Tests for module spectral_trace_list.py"""
import os
import pytest
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from synphot import Empirical1D, SourceSpectrum

from scopesim.effects import spectral_trace_list as spt
from scopesim.optics.fov_manager import FovVolumeList
from scopesim.effects.spectral_trace_list import SpectralTraceList, \
    SpectralTraceListWheel
from scopesim.effects.spectral_trace_list_utils import SpectralTrace
from scopesim.tests.mocks.py_objects import trace_list_objects as tlo
from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.base_classes import PoorMansHeader
from scopesim import rc

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SPEC/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]

if rc.__basic_inst_path__ not in rc.__search_path__:
    rc.__search_path__.append(rc.__basic_inst_path__)

PLOTS = False

# pylint: disable=missing-class-docstring,
# pylint: disable=missing-function-docstring


@pytest.fixture(name="slit_header", scope="class")
def fixture_slit_header():
    return ho._short_micado_slit_header()


@pytest.fixture(name="long_slit_header", scope="class")
def fixture_long_slit_header():
    return ho._long_micado_slit_header()


@pytest.fixture(name="full_trace_list", scope="class")
def fixture_full_trace_list():
    """Instantiate a trace definition hdu list"""
    return tlo.make_trace_hdulist()


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(SpectralTraceList(), SpectralTraceList)

    # @pytest.mark.usefixtures("full_trace_list")
    def test_initialises_with_a_hdulist(self, full_trace_list):
        spt = SpectralTraceList(hdulist=full_trace_list)
        assert isinstance(spt, SpectralTraceList)
        assert spt.get_data(2, fits.BinTableHDU)
        # next assert that dispersion axis determined correctly
        assert list(spt.spectral_traces.values())[2].dispersion_axis == 'y'

    def test_initialises_with_filename(self, mock_dir):
        micado_spec_dir = mock_dir / "MICADO_SPEC"
        with patch("scopesim.rc.__search_path__", [micado_spec_dir]):
            spt = SpectralTraceList(filename="TRACE_MICADO.fits",
                                    wave_colname="wavelength", s_colname="xi")
        assert isinstance(spt, SpectralTraceList)
        # assert that dispersion axis taken correctly from header keyword
        assert list(spt.spectral_traces.values())[2].dispersion_axis == 'y'

    def test_getitem_returns_spectral_trace(self, full_trace_list):
        slist = SpectralTraceList(hdulist=full_trace_list)
        assert isinstance(slist['Sheared'], SpectralTrace)

    def test_setitem_appends_correctly(self, full_trace_list):
        slist = SpectralTraceList(hdulist=full_trace_list)
        n_trace = len(slist.spectral_traces)
        spt = tlo.trace_1()
        slist["New trace"] = spt
        assert len(slist.spectral_traces) == n_trace + 1


@pytest.mark.skip(reason="Ignoring old Spectroscopy integration tests")
class TestGetFOVHeaders:
    @pytest.mark.usefixtures("full_trace_list", "slit_header")
    def test_gets_the_headers(self, full_trace_list, slit_header):
        sptl = spt.SpectralTraceList(hdulist=full_trace_list)
        params = {"pixel_scale": 0.015, "plate_scale": 0.26666,
                  "wave_min": 0.7, "wave_max": 2.5}
        hdrs = sptl.get_fov_headers(slit_header, **params)

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
        sptl = spt.SpectralTraceList(filename="TRACE_15arcsec.fits",
                                     s_colname="xi",
                                     wave_colname="lam",
                                     spline_order=1)
        params = {"wave_min": wave_min, "wave_max": wave_max,
                  "pixel_scale": 0.004, "plate_scale": 0.266666667}
        hdrs = sptl.get_fov_headers(slit_hdr, **params)
        assert isinstance(sptl, spt.SpectralTraceList)

        print(len(hdrs))

        if PLOTS:
            sptl.plot(wave_min, wave_max)

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
#         sptl = spt.SpectralTraceList(filename=fname, s_colname='xi')
#
#         fvl = FovVolumeList()
#         fvl = spt.apply_to(fvl)
#
#         assert len(fvl) == 17


@pytest.fixture(name="spectral_trace_list", scope="class")
def fixture_spectral_trace_list():
    """Instantiate a SpectralTraceList"""
    return SpectralTraceList(hdulist=tlo.make_trace_hdulist())


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


class TestRectification:
    def test_rectify_cube_not_implemented(self, spectral_trace_list):
        hdulist = fits.HDUList()
        with pytest.raises(NotImplementedError):
            spectral_trace_list.rectify_cube(hdulist)

    # def test_rectify_traces_needs_ximin_and_ximax(self, spectral_trace_list):
    #    hdulist = fits.HDUList([fits.PrimaryHDU()])
    #    with pytest.raises(KeyError):
    #        spectral_trace_list.rectify_traces(hdulist)


class TestSpectralTraceListWheel:
    @pytest.mark.usefixtures("no_file_error")
    def test_basic_init(self):
        """
        This is a super basic test just to see the thing basically works and
        parameters are passed correctly. Please feel free to improve this!!
        """
        kwargs = {"current_trace_list": "bogus",
                  "filename_format": "bogus_{}",
                  "trace_list_names": ["foo"]}
        stw = SpectralTraceListWheel(**kwargs)
        assert isinstance(stw, SpectralTraceListWheel)
        assert stw.meta["current_trace_list"] == "bogus"
        assert stw.meta["filename_format"] == "bogus_{}"
        assert stw.meta["trace_list_names"] == ["foo"]
        assert isinstance(stw.trace_lists["foo"], SpectralTraceList)
        assert stw.trace_lists["foo"].meta["filename"] == "bogus_foo"


class TestUnresolvedSpectralTraceListInit:
    def test_init(self):
        sptl = spt.UnresolvedSpectralTraceList(filename="INS_mos_traces.fits")
        assert len(sptl.spectral_traces) == 5
        assert sptl.spectral_traces["TRACE_Ap0"].meta["trace_id"] == "TRACE_Ap0"

        if PLOTS:
            for trace in list(sptl.spectral_traces.values()):
                plt.plot(trace.table["x"], trace.table["y"])
            plt.show()

    def test_get_xyz_for_trace(self):
        sptl = spt.UnresolvedSpectralTraceList(filename="INS_mos_traces.fits",
                                               plate_scale=20,
                                               pixel_scale=0.2)
        waves = np.linspace(1.75, 2.5, 101) * u.um
        flux = np.zeros_like(waves.value)
        flux[::10] = 1
        spec = SourceSpectrum(Empirical1D, points=waves, lookup_table=flux)

        for trace in sptl.spectral_traces.values():
            x, y, z = sptl.get_xyz_for_trace(trace, spec)
            plt.plot(x, y)

        if PLOTS:
            plt.show()
            plt.pause(0)

    def test_trace_to_image(self):
        sptl = spt.UnresolvedSpectralTraceList(filename="INS_mos_traces.fits",
                                               plate_scale=20,
                                               pixel_scale=0.2)
        waves = np.linspace(1.75, 2.5, 101) * u.um
        flux = np.ones_like(waves.value)
        flux[::10] = 10
        spec = SourceSpectrum(Empirical1D, points=waves, lookup_table=flux)

        for i, trace in enumerate(sptl.spectral_traces.values()):
            x, y, z = sptl.get_xyz_for_trace(trace, spec)
            image = sptl.trace_to_image(x, y, z)
            plt.subplot(2, 3, i+1)
            plt.imshow(image, origin="lower")

        if PLOTS:
            plt.pause(0)
            plt.show()

class TestMosaicSpectralTraceList():
    def test_intit(self):
        sptl = spt.MosaicSpectralTraceList(
            plate_scale=6.33333,
            pixel_scale=0.095,
            wave_min=1.420,
            wave_max=1.857,
        )
        assert len(sptl.spectral_traces) == 14
        print(sptl.spectral_traces["Trace_Ap0"].table)

        # assert sptl.spectral_traces["TRACE_Ap0"].meta["trace_id"] == "TRACE_Ap0"
