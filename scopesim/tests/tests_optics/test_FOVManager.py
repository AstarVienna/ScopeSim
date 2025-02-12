import pytest
from pytest import approx
import numpy as np

from scopesim.optics.fov_manager import FOVManager
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
        fov_volume = fovs[0].volume()

        assert len(fovs) == 1
        assert fov_volume["xs"][0] == approx(-1024 / 3600)      # [deg] 2k detector / pixel_scale
        assert fov_volume["waves"][0] == 0.6            # [um] filter blue edge

    @pytest.mark.parametrize("chunk_size, n_fovs",
                             [(500, 25), (512, 16), (1000, 9), (1024, 4)])
    def test_returns_n_fovs_for_smaller_chunk_size(self, chunk_size, n_fovs):
        effects = eo._mvs_effects_list()
        fov_man = FOVManager(effects=effects, pixel_scale=1, plate_scale=1,
                             max_segment_size=1024**2, chunk_size=1024)
        fovs = list(fov_man.generate_fovs_list())
        fov_volume = fovs[0].volume()

        assert len(fovs) == 4
        assert fov_volume["xs"][0] == approx(-1024 / 3600)     # [deg] 2k detector / pixel_scale
        assert fov_volume["waves"][0] == 0.6            # [um] filter blue edge

    def test_fov_volumes_have_detector_dimensions_from_detector_list(self):
        effects = eo._mvs_effects_list()
        fov_man = FOVManager(effects=effects, pixel_scale=1, plate_scale=1)
        _ = list(fov_man.generate_fovs_list())
        detector_limits = fov_man.volumes_list.detector_limits

        assert detector_limits["xd_min"] != 0.0
        assert detector_limits["yd_max"] != 0.0
