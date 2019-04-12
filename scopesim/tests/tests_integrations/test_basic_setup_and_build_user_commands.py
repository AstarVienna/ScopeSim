import pytest
import shutil
import os
import scopesim as sim
import scopesim.server.database as sim_db


_parent_path = "./test_downloads_dir/"


@pytest.fixture(scope="module")
def temp_directory_structure():
    # setup
    sim_db.rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"] = _parent_path
    sim_db.set_up_local_package_directory(_parent_path, True)

    sim_db.download_package("MICADO")
    sim_db.download_package("ELT")
    sim_db.download_package("MAORY_MCAO_4mas")

    # run tests
    yield

    # teardown
    shutil.rmtree(_parent_path)


@pytest.fixture(scope="class")
def cmd_micado():
    cmd = sim.UserCommands(instrument="MICADO", mode="MODE_MCAO_WIDE")
    return cmd


@pytest.mark.usefixtures("temp_directory_structure", "cmd_micado")
class TestWorkFlowForBuildingUserCommands:

    def test_make_user_commands_from_scratch(self, cmd_micado):
        assert isinstance(cmd_micado, sim.UserCommands)

    def test_psf_file_path_exists(self, cmd_micado):
        assert os.path.exists(cmd_micado["SCOPE_PSF_FILE"])

    def test_filter_file_path_exists(self, cmd_micado):
        cmd_micado.select_filter("J")
        assert os.path.exists(cmd_micado["INST_FILTER_TC"])

