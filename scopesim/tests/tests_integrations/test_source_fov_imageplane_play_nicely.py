# pylint: disable=no-self-use
# pylint: disable=missing-function-docstring
# pylint: disable=invalid-name

import numpy as np
import pytest
from pytest import approx

from astropy import units as u

from astropy.wcs import WCS

from scopesim.optics import image_plane_utils as imp_utils
from scopesim.optics.fov import FieldOfView2D
from scopesim.optics.fov_manager import chunk_edges
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


class TestInteractionBetweenSourceFOVImagePlane:
    """
    Test:
    - fov extracts correctly from source object and converts to image
    - fov image is correctly added to the image plane
    """

    def test_can_extract_the_source_in_a_fov(self, fov_hdr, comb_src,
                                             implane_hdr):

        fov = FieldOfView2D(fov_hdr, waverange=[0.5, 2.5]*u.um, area=1*u.m**2)
        imp = ImagePlane(implane_hdr)

        fov.extract_from(comb_src)
        fov.view()
        imp.add(fov.hdu, wcs_suffix="D")

        assert np.sum(imp.image) > 0

        if PLOTS:
            plt.subplot(131)
            plt.imshow(comb_src.fields[3].data, origin="lower", norm=LogNorm())

            plt.subplot(132)
            plt.imshow(fov.hdu.data, origin="lower", norm=LogNorm())

            plt.subplot(133)
            plt.imshow(imp.hdu.data, origin="lower", norm=LogNorm())
            plt.show()

def _checkerboard_source(naxis, pixel_scale, crpix=None):
    """0/1 checkerboard image source, one board square per pixel."""
    from astropy.io import fits
    from scopesim.source.source import Source
    from scopesim.source import source_templates as st

    nx, ny = naxis
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["LINEAR", "LINEAR"]
    wcs.wcs.cunit = ["deg", "deg"]
    wcs.wcs.cdelt = 2 * [pixel_scale / 3600]
    wcs.wcs.crval = [0, 0]
    wcs.wcs.crpix = crpix or [(nx + 1) / 2, (ny + 1) / 2]

    hdr = wcs.to_header()
    hdr["BUNIT"] = ""
    board = (np.indices((ny, nx)).sum(axis=0) % 2).astype(float)
    hdu = fits.ImageHDU(header=hdr, data=board)
    hdu.header["SPEC_REF"] = 0
    hdu.header["SPECMAG"] = 0
    hdu.header["SPECUNIT"] = "mag"
    return Source(image_hdu=hdu, spectra=[st.vega_spectrum()]), board


def _chunked_fov_headers(x_mm, y_mm, pixel_size, pixel_scale, chunk):
    """Reproduce FOVManager's chunking of a detector footprint."""
    plate_scale = pixel_scale / pixel_size
    lims = [v / pixel_size * pixel_scale
            for v in (x_mm[0], x_mm[1], y_mm[0], y_mm[1])]

    def edges(vmin, vmax):
        # the real chunk-edge helper, not a copy of it
        return [vmin, *chunk_edges(vmin, vmax, chunk * pixel_scale), vmax]

    xe, ye = edges(lims[0], lims[1]), edges(lims[2], lims[3])
    for x0, x1 in zip(xe[:-1], xe[1:]):
        for y0, y1 in zip(ye[:-1], ye[1:]):
            skyhdr = imp_utils.header_from_list_of_xy(
                [x0 / 3600, x1 / 3600], [y0 / 3600, y1 / 3600],
                pixel_scale=pixel_scale / 3600)
            dethdr, _ = imp_utils.det_wcs_from_sky_wcs(
                WCS(skyhdr), pixel_scale, plate_scale)
            skyhdr.update(dethdr.to_header())
            yield skyhdr


class TestChunkedFOVsTileTheImagePlane:
    """A checkerboard must survive being cut into chunked FOVs and reassembled.

    ``chunk_size`` is 8 here so the whole thing runs in well under a second,
    but the geometry is MICADO's: the detector footprint is not a whole number
    of pixels across (MICADO's y axis is 189.32 mm / 0.015 mm = 12621.33 px),
    so the last chunk in each axis is a fraction of a chunk wide and lands on
    a pixel lattice offset from its neighbours'. The projection used to snap
    those chunks inconsistently, which put a duplicated or a dropped
    row/column at the seam and flipped the checkerboard phase across it.
    """

    PIXEL_SIZE = 0.01     # mm / pix
    PIXEL_SCALE = 0.2     # arcsec / pix
    CHUNK = 8

    @pytest.mark.parametrize("half_mm, exact_extent", [
        (0.12,   24.0),    # whole pixels, whole chunks
        (0.105,  21.0),    # whole pixels, fractional last chunk
        (0.0866, 17.32),   # fractional pixels -- MICADO's y axis
        (0.0855, 17.10),
    ])
    def test_checkerboard_survives_chunking(self, half_mm, exact_extent):
        x_mm = y_mm = (-half_mm, half_mm)
        assert (x_mm[1] - x_mm[0]) / self.PIXEL_SIZE == approx(exact_extent)

        imphdr = imp_utils.header_from_list_of_xy(
            list(x_mm), list(y_mm), pixel_scale=self.PIXEL_SIZE, wcs_suffix="D")
        imp = ImagePlane(imphdr)

        source, board = _checkerboard_source(
            (imphdr["NAXIS1"], imphdr["NAXIS2"]), self.PIXEL_SCALE,
            crpix=[imphdr["CRPIX1D"], imphdr["CRPIX2D"]])

        fov_hdrs = list(_chunked_fov_headers(
            x_mm, y_mm, self.PIXEL_SIZE, self.PIXEL_SCALE, self.CHUNK))
        assert len(fov_hdrs) == 9

        for hdr in fov_hdrs:
            fov = FieldOfView2D(hdr, waverange=[0.5, 2.5] * u.um,
                                area=1 * u.m**2)
            fov.extract_from(source)
            fov.view()
            imp.add(fov.hdu, wcs_suffix="D")

        # no gaps and no double-counting at the chunk seams
        got = imp.data / imp.data.max()
        np.testing.assert_allclose(got, board, atol=1e-6)
