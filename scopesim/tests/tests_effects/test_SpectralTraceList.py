"""Tests for module spectral_trace_list.py"""

import pytest
from unittest.mock import patch

from astropy.io import fits


from scopesim.effects.spectral_trace_list import SpectralTraceList, \
    SpectralTraceListWheel
from scopesim.effects.spectral_trace_list_utils import SpectralTrace
from scopesim.tests.mocks.py_objects import trace_list_objects as tlo
from scopesim.tests.mocks.py_objects import header_objects as ho


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

    def test_initialises_with_a_hdulist(self, full_trace_list):
        spt = SpectralTraceList(hdulist=full_trace_list)
        assert isinstance(spt, SpectralTraceList)
        assert isinstance(spt.data_container._file[2], fits.BinTableHDU)
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


@pytest.fixture(name="spectral_trace_list", scope="class")
def fixture_spectral_trace_list():
    """Instantiate a SpectralTraceList"""
    return SpectralTraceList(hdulist=tlo.make_trace_hdulist())

#class TestRectification:

    # def test_rectify_traces_needs_ximin_and_ximax(self, spectral_trace_list):
    #    hdulist = fits.HDUList([fits.PrimaryHDU()])
    #    with pytest.raises(KeyError):
    #        spectral_trace_list.rectify_traces(hdulist)


class TestRectifyTracesInput:
    def test_rectify_traces_opens_filename_before_use(
            self, spectral_trace_list, tmp_path, monkeypatch):
        """rectify_traces accepts a filename, but passed the raw string on
        to make_image_interpolations and SpectralTrace.rectify instead of
        the opened HDUList."""
        import numpy as np
        from astropy.table import Table
        from scopesim.effects import spectral_trace_list as stl_mod

        readout = tmp_path / "readout.fits"
        fits.HDUList([fits.PrimaryHDU(),
                      fits.ImageHDU(data=np.zeros((10, 10)))]
                     ).writeto(readout)

        # Stub the filter lookup; it needs a fully configured instrument
        class FakeFilterCurve:
            table = Table(data=[[1.0, 2.0], [1.0, 1.0]],
                          names=["wavelength", "transmission"])

        monkeypatch.setattr(stl_mod, "FilterCurve",
                            lambda **kwargs: FakeFilterCurve())
        monkeypatch.setattr(stl_mod, "from_currsys",
                            lambda *args, **kwargs: "J")

        # Record what object each trace's rectify() receives
        received = []

        def fake_rectify(self, hdulist, **kwargs):
            received.append(hdulist)
            return None

        monkeypatch.setattr(SpectralTrace, "rectify", fake_rectify)

        result = spectral_trace_list.rectify_traces(
            str(readout), xi_min=-1, xi_max=1)

        assert isinstance(result, fits.HDUList)
        assert received, "rectify() was never called"
        assert all(isinstance(hdul, fits.HDUList) for hdul in received), \
            "traces received the filename string, not the opened HDUList"


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
