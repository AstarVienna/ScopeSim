
from copy import deepcopy
import pytest
from pytest import approx
from unittest.mock import patch

import numpy as np
from astropy import units as u
from astropy.table import Table

import scopesim as sim
from scopesim.optics.fov_manager import FOVManager
from scopesim.optics.image_plane import ImagePlane
from scopesim.optics.optical_train import OpticalTrain
from scopesim.optics.optics_manager import OpticsManager
from scopesim.optics.optical_element import OpticalElement
from scopesim.commands.user_commands import UserCommands
from scopesim.effects import Effect, DetectorList
from scopesim.utils import find_file

from scopesim.tests.mocks.py_objects import source_objects as src_objs

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


PLOTS = False


@pytest.fixture(scope="module", autouse=True)
def patch_globals():
    """Prevent modification of globals from within this module."""
    with patch("scopesim.rc.__currsys__"):
        yield


# TODO: check if class scope breaks anything (used to be function scope)
@pytest.fixture(scope="class")
def cmds(mock_path, mock_path_yamls):
    with patch("scopesim.rc.__search_path__", [mock_path, mock_path_yamls]):
        return UserCommands(yamls=[find_file("CMD_mvs_cmds.yaml")])

@pytest.fixture(scope="class")
def cmds_with_ignore(mock_path, mock_path_yamls):
    with patch("scopesim.rc.__search_path__", [mock_path, mock_path_yamls]):
        cmds = UserCommands(yamls=[find_file("CMD_mvs_cmds.yaml")])
        cmds.ignore_effects += ["detector QE curve"]
        return cmds


# TODO: check if class scope breaks anything (used to be function scope)
@pytest.fixture(scope="class")
def unity_cmds(mock_path, mock_path_yamls):
    with patch("scopesim.rc.__search_path__", [mock_path, mock_path_yamls]):
        return UserCommands(yamls=[find_file("CMD_unity_cmds.yaml")])


@pytest.fixture(scope="function")
def tbl_src():
    return src_objs._table_source()


@pytest.fixture(scope="function")
def im_src():
    return src_objs._image_source()


@pytest.fixture(scope="function")
def unity_src():
    return src_objs._unity_source(n=10001)


@pytest.fixture(scope="class")
def simplecado_opt(mock_path_yamls):
    simplecado_yaml = str(mock_path_yamls / "SimpleCADO.yaml")
    cmd = sim.UserCommands(yamls=[simplecado_yaml])
    return sim.OpticalTrain(cmd)


@pytest.mark.usefixtures("patch_mock_path")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(OpticalTrain(), OpticalTrain)

    def test_initialises_with_basic_commands(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt, OpticalTrain)

    def test_has_user_commands_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.cmds, UserCommands)

    def test_has_optics_manager_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.optics_manager, OpticsManager)

    def test_has_fov_manager_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        print(opt.fov_manager)
        assert isinstance(opt.fov_manager, FOVManager)

    def test_has_image_plane_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.image_planes[0], ImagePlane)

    def test_has_yaml_dict_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.yaml_dicts, list) and len(opt.yaml_dicts) > 0

    def test_ignore_effects_works(self, cmds_with_ignore):
        opt = OpticalTrain(cmds=cmds_with_ignore)
        assert opt["detector QE curve"].include is False


@pytest.mark.slow
@pytest.mark.usefixtures("patch_mock_path")
class TestObserve:
    """
    Almost all tests here are for visual inspection.
    No asserts, this just to test that the puzzle gets put back together
    after it is chopped up by the FOVs.
    """
    def test_observe_works_for_none(self, cmds):
        opt = OpticalTrain(cmds)
        opt.observe()
        empty = sim.source.source_templates.empty_sky()
        assert(opt._last_source.fields[0].field == empty.fields[0].field)

    def test_observe_works_for_table(self, cmds, tbl_src):
        opt = OpticalTrain(cmds)
        opt.observe(tbl_src)

        if PLOTS:
            plt.imshow(opt.image_planes[0].image.T, origin="lower",
                       norm=LogNorm())
            plt.show()

    def test_observe_works_for_image(self, cmds, im_src):
        opt = OpticalTrain(cmds)
        opt.observe(im_src)

        if PLOTS:
            plt.imshow(opt.image_planes[0].image.T, origin="lower",
                       norm=LogNorm())
            plt.show()

    def test_observe_works_for_source_distributed_over_several_fovs(self, cmds,
                                                                    im_src):

        orig_sum = np.sum(im_src.fields[0].data)

        cmds["SIM_PIXEL_SCALE"] = 0.02
        opt = OpticalTrain(cmds)
        opt.observe(im_src)

        wave = np.arange(0.5, 2.51, 0.1)*u.um
        unit = u.Unit("ph s-1 m-2 um-1")
        implane = opt.image_planes[0]
        final_sum = np.sum(implane.image)
        print(orig_sum, final_sum)

        if PLOTS:
            for fov in opt.fov_manager.fovs:
                cnrs = fov.corners[1]
                plt.plot(cnrs[0], cnrs[1])

            plt.imshow(implane.image.T, origin="lower", norm=LogNorm(),
                       extent=(-implane.hdu.header["NAXIS1"] / 2,
                               implane.hdu.header["NAXIS1"] / 2,
                               -implane.hdu.header["NAXIS2"] / 2,
                               implane.hdu.header["NAXIS2"] / 2,))
            plt.colorbar()
            plt.show()

    def test_observe_works_for_many_sources_distributed(self, cmds, im_src):
        orig_sum = np.sum(im_src.fields[0].data)
        im_src.fields[0].data += 1
        im_src1 = deepcopy(im_src)
        im_src2 = deepcopy(im_src)
        im_src2.shift(7, 7)
        im_src3 = deepcopy(im_src)
        im_src3.shift(-10, 14)
        im_src4 = deepcopy(im_src)
        im_src4.shift(-4, -6)
        im_src5 = deepcopy(im_src)
        im_src5.shift(15, -15)
        multi_img = im_src1 + im_src2 + im_src3 + im_src4 + im_src5

        cmds["SIM_PIXEL_SCALE"] = 0.02
        opt = OpticalTrain(cmds)
        opt.observe(multi_img)

        implane = opt.image_planes[0]
        final_sum = np.sum(implane.image)
        print(orig_sum, final_sum)

        if PLOTS:
            for fov in opt.fov_manager.fovs:
                cnrs = fov.corners[1]
                plt.plot(cnrs[0], cnrs[1])

            plt.imshow(implane.image.T, origin="lower", norm=LogNorm(),
                       extent=(-implane.hdu.header["NAXIS1"] / 2,
                               implane.hdu.header["NAXIS1"] / 2,
                               -implane.hdu.header["NAXIS2"] / 2,
                               implane.hdu.header["NAXIS2"] / 2,))
            plt.colorbar()
            plt.show()

    def test_works_with_a_pointer_to_fits_imagehdu(self, cmds):
        # Basically just checking to make sure observe doesn't throw an error
        # when passed a Source object with a file pointer ImageHDU
        fits_src = src_objs._fits_image_source()
        array_src = src_objs._image_source()

        src = fits_src + array_src
        opt = OpticalTrain(cmds)
        opt.observe(src)

        assert np.sum(opt.image_planes[0].data) > 0


@pytest.mark.usefixtures("patch_mock_path")
class TestReadout:
    def test_readout_works_when_source_observed(self, unity_cmds, unity_src):

        opt = OpticalTrain(unity_cmds)
        opt.observe(unity_src)
        hdus = opt.readout()
        hdu = hdus[0]

        if PLOTS:
            plt.subplot(221)
            plt.imshow(unity_src.fields[0].data)
            plt.colorbar()

            plt.subplot(222)
            plt.imshow(opt.image_planes[0].image)
            plt.colorbar()

            plt.subplot(223)
            plt.imshow(hdu[1].data)
            plt.colorbar()
            plt.show()

        _ = np.average(unity_src.fields[0].data)
        assert np.median(hdu[1].data) == approx(np.pi / 4., rel=1e-2)


class TestGetItems:
    def test_optical_element_returned_for_unique_name(self, simplecado_opt):
        print(type(simplecado_opt["test_detector"]))
        assert isinstance(simplecado_opt["test_detector"], OpticalElement)

    def test_effect_returned_for_unique_name(self, simplecado_opt):
        assert isinstance(simplecado_opt["test_detector_list"], Effect)

    def test_raises_error_for_bogus_string(self, simplecado_opt):
        with pytest.raises(ValueError):
            simplecado_opt["bogus"]

    def test_list_of_effects_returned_for_effect_class(self, simplecado_opt):
        effects = simplecado_opt[DetectorList]
        assert isinstance(effects, list)
        assert len(effects) == 1


class TestSetItems:
    def test_effect_kwarg_can_be_changed_using(self, simplecado_opt):
        simplecado_opt["dark_current"] = {"value": 0.2}
        assert simplecado_opt["dark_current"].meta["value"] == 0.2

    def test_effect_include_can_be_toggled_with_setitem(self, simplecado_opt):
        assert simplecado_opt["dark_current"].include is True
        simplecado_opt["dark_current"].include = False
        assert simplecado_opt["dark_current"].include is False
        assert simplecado_opt["dark_current"].meta["include"] is False


class TestListEffects:
    def test_effects_listed_in_table(self, simplecado_opt):
        assert isinstance(simplecado_opt.effects, Table)
        simplecado_opt["dark_current"].include = False
        assert bool(simplecado_opt.effects["included"][1]) is False
        simplecado_opt["alt_dark_current"].include = True
        assert bool(simplecado_opt.effects["included"][2]) is True

        print("\n", simplecado_opt.effects)


class TestShutdown:
    """Test that fits files are closed on shutdown of OpticalTrain"""

    def test_files_closed_on_shutdown(self, simplecado_opt, mock_path):
        """Test for closed files in two ways:
        - `closed` flag is set to True
        - data access fails
        """
        # Add an effect with a psf
        with patch("scopesim.rc.__search_path__", [mock_path]):
            psf = sim.effects.FieldConstantPSF(filename="test_ConstPSF.fits",
                                               name="testpsf",
                                               cmds=simplecado_opt.cmds)
        simplecado_opt.optics_manager.add_effect(psf)
        # This is just to make sure that we have an open file
        assert not simplecado_opt['testpsf']._file._file.closed

        simplecado_opt.shutdown()

        # 1. Check the `closed` flags where available
        flags = []
        for effect_name in simplecado_opt.effects['name']:
            try:
                flags.append(simplecado_opt[effect_name]._file._file.closed)
            except AttributeError:
                pass

        assert all(flags)

        # 2. Check that data access fails
        with pytest.raises(ValueError):
            print(simplecado_opt['testpsf']._file[2].data)
