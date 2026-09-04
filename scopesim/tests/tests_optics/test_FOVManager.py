import pytest
from pytest import approx
import numpy as np
from astropy import units as u

from scopesim.optics.fov_manager import FOVManager, chunk_edges
from scopesim.tests.mocks.py_objects import effects_objects as eo
from scopesim.utils import from_currsys


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(preload_fovs=False), FOVManager)

    @pytest.mark.usefixtures("patch_mock_path")
    def test_initialises_with_list_of_effects(self):
        effects = eo._mvs_effects_list()
        assert isinstance(FOVManager(effects, preload_fovs=False), FOVManager)


@pytest.mark.usefixtures("patch_mock_path")
class TestGenerateFovList:
    @pytest.mark.slow
    def test_returns_default_single_entry_fov_list_for_no_effects(self):
        fov_man = FOVManager(pixel_scale=1, plate_scale=1)
        assert len(fov_man.volumes_list) == 1, "volumes_list should have only 1 element initially."
        fov_vol_org = fov_man.volumes_list[0]

        chunk_size = from_currsys(fov_man.meta["chunk_size"], fov_man.cmds)
        n_vol_x = len(np.arange(fov_vol_org["x_min"], fov_vol_org["x_max"], chunk_size))
        n_vol_y = len(np.arange(fov_vol_org["y_min"], fov_vol_org["y_max"], chunk_size))
        fovs = list(fov_man.generate_fovs_list())

        assert len(fovs) == n_vol_x * n_vol_y, (f"Expected {n_vol_x} * {n_vol_y} = {n_vol_x * n_vol_y} volumes, "
                                                f"but got {len(fovs)} volumes")

    def test_returns_single_fov_for_mvs_system(self):
        effects = eo._mvs_effects_list()
        fov_man = FOVManager(effects=effects, pixel_scale=1, plate_scale=1)
        fovs = list(fov_man.generate_fovs_list())
        fov_skycorners, _ = fovs[0].get_corners("deg")

        assert len(fovs) == 1
        assert fov_skycorners.min(axis=0)[0] == approx(-1024 / 3600)  # [deg] 2k detector / pixel_scale
        assert fovs[0].waverange[0] == 0.6 * u.um  # filter blue edge

    @pytest.mark.parametrize("chunk_size, n_fovs",
                             [(500, 25), (512, 16), (1000, 9), (1024, 4)])
    def test_returns_n_fovs_for_smaller_chunk_size(self, chunk_size, n_fovs):
        effects = eo._mvs_effects_list()
        fov_man = FOVManager(effects=effects, pixel_scale=1, plate_scale=1,
                             max_segment_size=1024**2, chunk_size=1024)
        fovs = list(fov_man.generate_fovs_list())
        fov_skycorners, _ = fovs[0].get_corners("deg")

        assert len(fovs) == 4
        assert fov_skycorners.min(axis=0)[0] == approx(-1024 / 3600)  # [deg] 2k detector / pixel_scale
        assert fovs[0].waverange[0] == 0.6 * u.um  # filter blue edge


class TestChunkEdges:
    """Chunk edges must be interior and platform-independent.

    np.arange(vmin, vmax, step) returns one element too many when
    (vmax - vmin) / step rounds just above an integer, and whether that
    element lands on vmax or one ulp below it depends on how the platform
    rounds vmin + i * step. One ulp below vmax is a legal split point, so
    FovVolumeList.split would carve off a ~1e-15 px volume that
    create_wcs_from_points then rounds up into a spurious one-pixel FOV.
    """

    @pytest.mark.parametrize("span_px, chunk, n_chunks", [
        (24.0, 8, 3),          # exact multiple -- the fragile case
        (24.0, 4, 6),          # exact multiple
        (16.0, 8, 2),
        (21.0, 8, 3),          # fractional last chunk
        (17.32, 8, 3),         # fractional pixels, fractional last chunk
        (12621 + 1 / 3, 2048, 7),   # MICADO's y axis
        (5.0, 8, 1),           # smaller than one chunk
    ])
    def test_edges_are_interior_and_count_is_right(self, span_px, chunk,
                                                   n_chunks):
        pixel_scale = 0.1
        vmin = -span_px / 2 * pixel_scale
        vmax = vmin + span_px * pixel_scale
        edges = chunk_edges(vmin, vmax, chunk * pixel_scale)

        assert len(edges) == n_chunks - 1
        assert (edges > vmin).all()
        assert (edges < vmax).all()

    @pytest.mark.parametrize("span_px, chunk", [(24.0, 8), (24.0, 4),
                                                (16.0, 8), (32.0, 16)])
    def test_exact_multiples_give_no_sliver_chunk(self, span_px, chunk):
        # A span that is an exact multiple of the chunk must come out as
        # equal chunks, with no degenerate remainder however the float
        # arithmetic rounds.
        pixel_scale = 0.1
        vmin = -span_px / 2 * pixel_scale
        vmax = vmin + span_px * pixel_scale
        edges = chunk_edges(vmin, vmax, chunk * pixel_scale)

        bounds = np.array([vmin, *edges, vmax])
        widths = np.diff(bounds) / pixel_scale
        assert len(widths) == span_px / chunk
        np.testing.assert_allclose(widths, chunk, atol=1e-9)

    def test_a_split_one_ulp_inside_the_upper_bound_is_not_emitted(self):
        # np.arange does emit such a value on arm64 macOS; assert directly
        # that chunk_edges never does, on any platform.
        vmin, vmax, step = -1.2000000000000002, 1.2000000000000002, 0.8
        edges = chunk_edges(vmin, vmax, step)
        assert (edges < np.nextafter(vmax, -np.inf)).all()
        # ... which np.arange does not guarantee
        assert len(np.arange(vmin, vmax, step)) == len(edges) + 2
