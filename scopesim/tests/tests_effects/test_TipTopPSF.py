# -*- coding: utf-8 -*-
"""Tests for the TipTopPSF effect.

The TipTop server is never contacted: the tiptop_ipy query function is
monkeypatched to return a synthetic Gaussian PSF.
"""

import numpy as np
import pytest
from astropy.io import fits

from scopesim.effects import TipTopPSF
from scopesim.effects.psfs import tiptop as tiptop_module


PIXEL_SCALE = 0.004     # [arcsec]
N_PIX = 128


class FakeResult:
    """Mimics tiptop_ipy.TipTopResult for a single on-axis PSF."""

    def __init__(self, n_pix):
        x = np.arange(n_pix) - n_pix / 2 + 0.5
        xx, yy = np.meshgrid(x, x)
        self.psf = np.exp(-(xx**2 + yy**2) / (2 * 3.0**2))[None, :, :]
        self.strehl = np.array([0.42])
        self.fwhm = np.array([12.3])


@pytest.fixture(name="mock_server")
def fixture_mock_server(monkeypatch):
    """Replace the server query with a local Gaussian PSF factory."""
    calls = []

    def fake_generate(self, conn, path, wavelength, n_pix, pixel_scale):
        calls.append(wavelength)
        result = FakeResult(n_pix)
        kernel = result.psf[0] / result.psf[0].sum()
        hdu = fits.ImageHDU(kernel.astype(np.float32))
        hdu.header["WAVE0"] = wavelength
        hdu.header["CDELT1"] = pixel_scale
        hdu.header["CDELT2"] = pixel_scale
        hdu.header["CUNIT1"] = "arcsec"
        hdu.header["CUNIT2"] = "arcsec"
        hdu.header["STREHL"] = float(result.strehl[0])
        path.parent.mkdir(parents=True, exist_ok=True)
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)

    monkeypatch.setattr(TipTopPSF, "_generate_psf_file", fake_generate)
    return calls


def make_psf(tmp_path, **kwargs):
    params = {
        "instrument": "MAVIS",
        "wavelength": 0.55,
        "fov_pix": N_PIX,
        "pixel_scale": PIXEL_SCALE,
        "cache_dir": str(tmp_path / "cache"),
    }
    params.update(kwargs)
    return TipTopPSF(**params)


class TestInit:
    def test_initialises_without_server_contact(self, tmp_path):
        psf = make_psf(tmp_path)
        assert isinstance(psf, TipTopPSF)
        assert not (tmp_path / "cache").exists()

    def test_throws_without_instrument_or_ini(self):
        with pytest.raises(ValueError):
            TipTopPSF(wavelength=0.55)

    def test_clips_fov_to_max(self, tmp_path, mock_server):
        psf = make_psf(tmp_path, fov_pix=4096)
        n_pix, _ = psf._grid_parameters()
        assert n_pix == tiptop_module.MAX_PSF_SIDE


class TestGetKernel:
    def test_returns_normalised_kernel(self, tmp_path, mock_server):
        psf = make_psf(tmp_path)
        kernel = psf.get_kernel(PIXEL_SCALE)
        assert kernel.sum() == pytest.approx(1, rel=1e-6)
        assert kernel.shape == (N_PIX, N_PIX)

    def test_generates_only_once_per_wavelength(self, tmp_path, mock_server):
        psf = make_psf(tmp_path)
        psf.get_kernel(PIXEL_SCALE)
        psf.get_kernel(PIXEL_SCALE)
        assert len(mock_server) == 1

    def test_second_instance_uses_disk_cache(self, tmp_path, mock_server):
        psf1 = make_psf(tmp_path)
        kernel1 = psf1.get_kernel(PIXEL_SCALE)
        psf2 = make_psf(tmp_path)
        kernel2 = psf2.get_kernel(PIXEL_SCALE)
        assert len(mock_server) == 1
        assert np.allclose(kernel1, kernel2)

    def test_rescales_to_fov_pixel_scale(self, tmp_path, mock_server):
        psf = make_psf(tmp_path)
        kernel = psf.get_kernel(2 * PIXEL_SCALE)
        assert kernel.shape[0] == pytest.approx(N_PIX / 2, abs=1)
        assert kernel.sum() == pytest.approx(1, rel=1e-6)

    def test_strehl_ratio_from_cache(self, tmp_path, mock_server):
        psf = make_psf(tmp_path)
        psf.get_kernel(PIXEL_SCALE)
        assert psf.strehl_ratio == pytest.approx(0.42)


class TestWaveDict:
    def test_filter_name_resolves_wavelength(self, tmp_path, mock_server):
        psf = make_psf(tmp_path, wavelength=None,
                       wave_dict={"V": 0.545, "R": 0.641},
                       filter_name="V")
        psf.get_kernel(PIXEL_SCALE)
        assert mock_server == [0.545]

    def test_unknown_filter_raises(self, tmp_path, mock_server):
        psf = make_psf(tmp_path, wavelength=None,
                       wave_dict={"V": 0.545}, filter_name="K")
        with pytest.raises(KeyError):
            psf.get_kernel(PIXEL_SCALE)

    def test_one_call_per_filter(self, tmp_path, mock_server):
        for filter_name in ("V", "R", "V", "R"):
            psf = make_psf(tmp_path, wavelength=None,
                           wave_dict={"V": 0.545, "R": 0.641},
                           filter_name=filter_name)
            psf.get_kernel(PIXEL_SCALE)
        assert sorted(mock_server) == [0.545, 0.641]
