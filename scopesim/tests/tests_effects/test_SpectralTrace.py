import numpy as np
import pytest
from _pytest.python_api import approx
from astropy.wcs import WCS
from matplotlib import pyplot as plt

from scopesim.effects.spectral_trace_list_utils import SpectralTrace
from scopesim.effects.spectral_trace_list_utils import get_affine_parameters
from scopesim.tests.mocks.py_objects import trace_list_objects as tlo
from scopesim.tests.tests_effects.test_SpectralTraceList import PLOTS

PLOTS = False



@pytest.fixture(scope="function")
def basic_trace():
    return tlo.trace_0()


@pytest.fixture(scope="function")
def horizontal_trace():
    return tlo.trace_4()


@pytest.fixture(scope="function")
def diagonal_trace():
    return tlo.trace_2()


@pytest.fixture(scope="function")
def curved_trace():
    return tlo.trace_3()


class TestInit:
    @pytest.mark.usefixtures("basic_trace")
    def test_initialised_with_curve(self, basic_trace):
        assert isinstance(SpectralTrace(basic_trace), SpectralTrace)

    def test_initialised_with_nothing_throws_error(self):
        with pytest.raises(TypeError):
            SpectralTrace()


class TestGetMaxDispersion:
    @pytest.mark.usefixtures("basic_trace")
    def test_dispersion_for_vertically_aligned_trace(self, basic_trace):
        # ..todo: accuracy of get_max_dispersion should be tested in a trace_list_utils tests file
        spt = SpectralTrace(basic_trace)
        disp, wave = spt.get_max_dispersion()
        # dispersion is calculated by distance [mm] / wavelength coverage [um]
        dy = np.diff(basic_trace["y"])
        dw = np.diff(basic_trace["wavelength"])
        assert np.average(disp) == approx(np.average(dy / dw), rel=1e-5)
        assert len(disp) == len(wave)
        assert all(np.diff(wave) > 0)

    @pytest.mark.usefixtures("horizontal_trace")
    def test_dispersion_for_horizontally_aligned_trace(self, horizontal_trace):
        spt = SpectralTrace(horizontal_trace)
        disp, wave = spt.get_max_dispersion()
        dx = np.abs(np.diff(horizontal_trace["x"]))
        dw = np.abs(np.diff(horizontal_trace["wavelength"]))

        assert np.average(disp) == approx(np.average(dx / dw), rel=1e-5)

    def test_dispersion_for_curved_traces(self):
        # .. todo: come back to this for curved traces that go beyond 45 deg
        pass


class TestGetPixelEdges:
    @pytest.mark.usefixtures("basic_trace")
    def test_monochromatic_trace_curves_are_one_pixel_apart(self, basic_trace):
        pixel_size = 0.015
        spt = SpectralTrace(basic_trace)
        disp, wave = spt.get_max_dispersion()
        wbedges = spt.get_pixel_wavelength_edges(pixel_size)     # wavelength edges
        dist_between_mtc = np.average(disp) * np.average(np.diff(wbedges))

        assert dist_between_mtc == approx(pixel_size, rel=1e-5)

    @pytest.mark.usefixtures("horizontal_trace")
    def test_mtc_one_pixel_apart_for_horizontal_traces(self, horizontal_trace):
        pixel_size = 0.015
        spt = SpectralTrace(horizontal_trace)
        disp, wave = spt.get_max_dispersion()
        wbedges = spt.get_pixel_wavelength_edges(pixel_size)  # wavelength edges
        dist_between_mtc = np.average(disp) * np.average(np.diff(wbedges))

        assert dist_between_mtc == approx(pixel_size, rel=1e-5)

    @pytest.mark.usefixtures("diagonal_trace")
    def test_mtc_one_pixel_apart_for_diagonal_traces(self, diagonal_trace):
        pixel_size = 0.015
        spt = SpectralTrace(diagonal_trace)
        disp, wave = spt.get_max_dispersion()
        wave_edges = spt.get_pixel_wavelength_edges(pixel_size)  # wavelength edges
        dist_between_mtc = np.average(disp) * np.average(np.diff(wave_edges))

        assert dist_between_mtc == approx(pixel_size, rel=1e-5)


@pytest.mark.skip(reason="MonochromaticTraceCurves not needed for SpecCADO")
class TestGetTraceCurves:
    @pytest.mark.usefixtures("basic_trace")
    def test_mtc_are_one_pixel_removed_from_each_other(self, basic_trace):
        spt = SpectralTrace(basic_trace)
        mtcs = spt.get_trace_curves(0.015)
        pix_cen_x = [mtc.header["CRVAL1D"]/0.015 for mtc in mtcs]
        pix_cen_y = [mtc.header["CRVAL2D"]/0.015 for mtc in mtcs]
        assert np.average(np.diff(pix_cen_y)) == approx(1, rel=1e-5)

    @pytest.mark.usefixtures("horizontal_trace")
    def test_mtc_distances_are_one_pixel_horiz_trace(self, horizontal_trace):
        spt = SpectralTrace(horizontal_trace)
        mtcs = spt.get_trace_curves(0.015)
        pix_cen_x = [mtc.header["CRVAL1D"]/0.015 for mtc in mtcs]
        pix_cen_y = [mtc.header["CRVAL2D"]/0.015 for mtc in mtcs]
        assert np.abs(np.average(np.diff(pix_cen_x))) == approx(1, rel=1e-5)

        if PLOTS:
            plt.plot(pix_cen_x, pix_cen_y)
            plt.show()

    @pytest.mark.usefixtures("diagonal_trace")
    def test_mtc_distances_are_one_pixel_diagonal_trace(self, diagonal_trace):
        # diagonal trace is 30 degrees off vertical
        spt = SpectralTrace(diagonal_trace)
        mtcs = spt.get_trace_curves(0.015)
        pix_cen_x = [mtc.header["CRVAL1D"] / 0.015 for mtc in mtcs]
        pix_cen_y = [mtc.header["CRVAL2D"] / 0.015 for mtc in mtcs]
        assert np.abs(np.average(np.diff(pix_cen_y))) == approx(1, rel=1e-5)

        if PLOTS:
            plt.plot(pix_cen_x, pix_cen_y, "o")
            plt.show()

    @pytest.mark.usefixtures("diagonal_trace")
    def test_limits_trace_curves_to_xy_edges(self, diagonal_trace):
        spt = SpectralTrace(diagonal_trace)
        xy_edges = {"x_min": -25, "x_max": -15, "y_min": 10, "y_max": 20}
        mtcs_all = spt.get_trace_curves(0.015)
        mtcs_xy_limited = spt.get_trace_curves(0.015, xy_edges=xy_edges)

        assert len(mtcs_all) > len(mtcs_xy_limited)

        if PLOTS:
            for mtc in mtcs_xy_limited:
                plt.plot(mtc.x, mtc.y)
            spt.plot(spt.wave_min, spt.wave_max)
            plt.axhline(xy_edges["y_min"])
            plt.axhline(xy_edges["y_max"])
            plt.axvline(xy_edges["x_min"])
            plt.axvline(xy_edges["x_max"])
            plt.show()


@pytest.mark.skip(reason="MonochromaticTraceCurves not needed for SpecCADO")
class TestGetCurveHeaders:
    @pytest.mark.usefixtures("basic_trace")
    def test_vertical_headers_are_all_one_pixel_apart(self, basic_trace):
        spt = SpectralTrace(basic_trace)
        pixel_size = 0.015
        hdrs = spt.get_curve_headers(pixel_size)
        x = [hdr["CRVAL1D"] / pixel_size for hdr in hdrs]
        y = [hdr["CRVAL2D"] / pixel_size for hdr in hdrs]
        assert np.all(np.abs(np.diff(y)) == approx(1, rel=1e-5))

    @pytest.mark.usefixtures("diagonal_trace")
    def test_diagonal_headers_are_all_one_pixel_apart(self, diagonal_trace):
        spt = SpectralTrace(diagonal_trace)
        pixel_size = 0.015
        hdrs = spt.get_curve_headers(pixel_size)
        x = [hdr["CRVAL1D"] / pixel_size for hdr in hdrs]
        y = [hdr["CRVAL2D"] / pixel_size for hdr in hdrs]
        assert np.all(np.abs(np.diff(y)) == approx(1, rel=1e-5))

    @pytest.mark.usefixtures("horizontal_trace")
    def test_diagonal_headers_are_all_one_pixel_apart(self, horizontal_trace):
        spt = SpectralTrace(horizontal_trace)
        pixel_size = 0.015
        hdrs = spt.get_curve_headers(pixel_size)
        x = [hdr["CRVAL1D"] / pixel_size for hdr in hdrs]
        y = [hdr["CRVAL2D"] / pixel_size for hdr in hdrs]
        assert np.all(np.abs(np.diff(x)) == approx(1, rel=1e-5))

    @pytest.mark.usefixtures("curved_trace")
    def test_curved_headers_are_all_one_pixel_apart(self, curved_trace):
        curved_trace["y1"] *= 1.1
        curved_trace["y2"] *= 1.2
        spt = SpectralTrace(curved_trace)
        pixel_size = 0.015
        hdrs = spt.get_curve_headers(pixel_size)
        dx = np.diff([hdr["CRVAL1D"] for hdr in hdrs])
        dy = np.diff([hdr["CRVAL2D"] for hdr in hdrs])
        dr = (dx**2 + dy**2)**0.5
        assert np.all(dr <= 1)

        # !!! PLOT this again to see issues
        if PLOTS:
            # orig world coords
            for row in curved_trace:
                x = [row["x0"], row["x2"]]
                y = [row["y0"], row["y2"]]
                len_mm = (np.diff(x)**2 + np.diff(y)**2)**0.5
                plt.plot(x, y, "k")
                plt.plot(x[0], y[0], "ko")
                plt.text(x[0], y[0], len_mm)

            # pixel coords
            for hdr in hdrs[::86]:
                xp = [0, hdr["NAXIS1"]]
                yp = [0, hdr["NAXIS2"]]
                wcs = WCS(hdr, key="D")
                # world coords
                xw, yw = wcs.all_pix2world(xp, yp, 1)
                plt.plot(xw, yw, "r")
                plt.plot(hdr["CRVAL1D"], hdr["CRVAL2D"], "ro")
                len_mm = (np.diff(xw)**2 + np.diff(yw)**2)**0.5
                plt.text(hdr["CRVAL1D"], hdr["CRVAL2D"], len_mm, color="red")
            plt.show()


class TestGetAffineParameters:
    @pytest.mark.parametrize("rot_ang, shear_ang",
                             [(0, 0), (45, 45), (60, 30), (-5, -15)]
                             )
    def test_returns_proper_shear_and_rotation_for_square(self, rot_ang,
                                                          shear_ang):
        x0 = -np.ones(5)
        x1 = np.ones(5)
        y0 = np.arange(5)
        y1 = np.arange(5)

        t = np.tan(shear_ang / 57.29578)
        x0 = x0 + t * y0
        x1 = x1 + t * y1

        c = np.cos(rot_ang / 57.29578)
        s = np.sin(rot_ang / 57.29578)
        x0a = c * x0 + -s * y0
        y0a = s * x0 + c * y0
        x1a = c * x1 + -s * y1
        y1a = s * x1 + c * y1

        coords = {"x": np.array([x0a, x1a]),
                  "y": np.array([y0a, y1a])}
        rot, shear = get_affine_parameters(coords)

        assert np.average(rot) == approx(rot_ang, abs=1e-5)
        assert -np.average(shear) == approx(shear_ang, abs=1e-5)
        # ..todo:: work out why this is negative? Definitionssache?

        if PLOTS:
            plt.figure(figsize=(6, 6))
            plt.plot(coords["x"][0], coords["y"][0])
            plt.plot(coords["x"][-1], coords["y"][-1])
            plt.xlim(-5, 5)
            plt.ylim(-5, 5)
            plt.show()

    @pytest.mark.skip(reason="no way of currently testing this")
    @pytest.mark.usefixtures("curved_trace")
    def test_no_same_angles_for_curved_trace(self, curved_trace):
        spt = SpectralTrace(curved_trace)
        mtcs = spt.get_trace_curves(0.015)
        rots = [mtc.meta["rotation"] for mtc in mtcs]
        shears = [mtc.meta["shear"] for mtc in mtcs]

        assert len(np.unique(rots)) == len(np.unique(np.diff(rots)))

        if PLOTS:
            plt.subplot(121)
            for mtc in mtcs[::100]:
                plt.plot(mtc.x, mtc.y)

            plt.subplot(122)
            plt.plot(rots)
            plt.plot(shears)

            plt.show()


class TestPlot:
    @pytest.mark.usefixtures("curved_trace")
    def test_plots(self, curved_trace):
        spt = SpectralTrace(curved_trace)
        if PLOTS:
            spt.plot(0.5, 2.5)
            plt.show()


class TestRepr:
    @pytest.mark.usefixtures("curved_trace")
    def test_list_info(self, curved_trace):
        spt = SpectralTrace(curved_trace, aperture_id=0, extension_id=2,
                            image_plane_id=0)
        print(spt)
