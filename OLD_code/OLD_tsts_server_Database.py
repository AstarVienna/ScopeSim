import os
import shutil
import pytest
from astropy.table import Table
import OLD_code.OLD_database as sim_db
import scopesim.utils

_parent_path = "./test_downloads_dir/"

# .. todo:: add test for what happens when there is no connection to the server


@pytest.fixture(scope="module")
def temp_directory_structure():
    # setup
    sim_db.rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"] = _parent_path
    sim_db.set_up_local_package_directory(_parent_path, True)

    # run tests
    yield

    # teardown
    shutil.rmtree(_parent_path)


@pytest.fixture(scope="class")
def download_test_package():
    return sim_db.download_package("test_package")


@pytest.mark.usefixtures("temp_directory_structure")
class TestGetLocalPackages:
    def test_returns_table_if_local_db_file_exists(self):
        local_dbs = sim_db.set_local_path_names(_parent_path)
        for path in local_dbs:
            db_tbl = sim_db.get_local_packages(path)
            assert isinstance(db_tbl, Table)

    def test_throws_error_if_local_db_file_path_is_bogus(self):
        local_db_path = "bogus.txt"
        with pytest.raises(ValueError):
            sim_db.get_local_packages(local_db_path)

    def test_empty_table_has_column_type_string(self):
        pass


class TestGetServerPackages:
    def test_returns_none_for_wrong_path(self):
        svr_db_url = "www.my-server.bogus"
        assert sim_db.get_server_packages(svr_db_url) is None

    def test_returns_table_if_path_correct(self):
        svr_path = sim_db._svr_inst_db()
        assert type(sim_db.get_server_packages(svr_path)) == Table
        svr_path = sim_db._svr_psf_db()
        assert type(sim_db.get_server_packages(svr_path)) == Table
        svr_path = sim_db._svr_src_db()
        assert type(sim_db.get_server_packages(svr_path)) == Table


class TestRenameServerTable:
    def test_throws_exception_when_not_enough_column_names(self):
        svr_table = Table(data=[[1, 1], [1, 1]])
        svr_table.meta["comments"] = ["# name"]

        with pytest.raises(Exception):
            sim_db.rename_table_colnames(svr_table)

    def test_renames_when_num_cols_equals_num_names(self):
        svr_table = Table(data=[[1, 1], [1, 1]])
        svr_table.meta["comments"] = ["# name other"]

        renamed_tbl = sim_db.rename_table_colnames(svr_table)
        assert renamed_tbl.colnames[0] == "name"
        assert renamed_tbl.colnames[1] == "other"


class TestCheckPackageExists:
    def test_returns_true_for_package_name(self):
        assert sim_db.check_package_exists("test_package") is True

    def test_returns_false_for_bogus_pkg_name(self):
        assert sim_db.check_package_exists("bogus") is False

    def test_returns_exception_if_package_path_is_broken(self):
        with pytest.raises(ValueError):
            sim_db.check_package_exists("non_existent_pkg")


class TestGetServerPackagePath:
    def test_return_url_for_existing_package(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        path = sim_db.get_server_package_path("test_package", svr_table)
        assert path == "test_package.zip"

    def test_returns_none_if_package_not_in_dbs(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        path = sim_db.get_server_package_path("bogus", svr_table)
        assert path is None


class TestGetPackageTableEntry:
    def test_returns_url_with_right_data(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        return_tbl = sim_db.get_package_table_entry("test_package", svr_table)
        assert return_tbl["path"] == "test_package.zip"

    def test_returns_none_if_package_not_it_table(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        assert sim_db.get_package_table_entry("bogus", svr_table) is None

    def test_returns_newest_with_multiple_entries(self):
        pass

    def test_returns_oldest_with_multiple_entries(self):
        pass


@pytest.mark.usefixtures("temp_directory_structure")
class TestDownloadPackage:
    def test_raise_error_when_pkg_not_in_db(self):
        with pytest.raises(ValueError):
            sim_db.download_package("bogus")

    def test_raise_error_when_pkg_file_doesnt_exist(self):
        with pytest.raises(ValueError):
            sim_db.download_package("non_existent_pkg")

    # avoid a test that is dependent on the network
    # ::todo add this to the integration test suite

    def test_package_file_exists_on_local_drive(self):
        local_filename = sim_db.download_package("test_package")
        assert os.path.exists(local_filename)

    def test_package_added_to_local_db(self):
        local_pkgs_before = sim_db.get_local_packages(sim_db._local_inst_db())
        sim_db.download_package("test_package")
        local_pkgs_after = sim_db.get_local_packages(sim_db._local_inst_db())
        assert len(local_pkgs_after) == len(local_pkgs_before) + 1


@pytest.mark.usefixtures("temp_directory_structure")
class TestAddPackageToLocalDb:
    def test_adds_row(self):
        local_table = Table(names=["name", "author", "date_added",
                                   "date_modified", "path"],
                            data=[["test_package"], ["Kieran Leschinski"],
                                  ["2018-11-09"], ["2018-11-09"],
                                  ["test_package.zip"]])
        len_local_table = len(local_table)
        pkg_entry = local_table[0]
        new_local_table = sim_db.add_pkg_to_local_db(pkg_entry, local_table)

        assert len(new_local_table) == len_local_table + 1

    def test_throws_exception_if_pkg_entry_isnt_astropy_row_class(self):
        local_table = Table(names=["name", "path"], data=[["test_package"],
                                                          ["test_package.zip"]])
        with pytest.raises(ValueError):
            sim_db.add_pkg_to_local_db("hello world!", local_table)

    def test_existing_is_renamed_if_new_package_is_newer(self):
        local_table = Table(names=["name", "date_modified"],
                            data=[["test_package"], ["2018-11-09"]])
        pkg_entry = Table(names=["name", "date_modified"],
                          data=[["test_package"], ["2018-11-14"]])
        new_local_table = sim_db.add_pkg_to_local_db(pkg_entry[0], local_table)
        assert new_local_table["name"][0] == "test_package_2018-11-09"
        assert new_local_table["name"][1] == "test_package"

    def test_new_package_is_renamed_if_new_package_is_older(self):
        local_table = Table(names=["name", "date_modified"],
                            data=[["test_package"], ["2018-11-14"]])
        pkg_entry = Table(names=["name", "date_modified"],
                          data=[["test_package"], ["2018-11-09"]])
        new_local_table = sim_db.add_pkg_to_local_db(pkg_entry[0], local_table)
        assert new_local_table["name"][0] == "test_package"
        assert new_local_table["name"][1] == "test_package_2018-11-09"


class TestChangeTableEntry:
    def test_string_changed_successfully_for_pattern(self):
        tbl = Table(names=["id", "name"], data=[[0], ["seb skelly"]])
        tbl = scopesim.utils.change_table_entry(tbl, "name", "seb skelly rocks",
                                        "seb skelly")
        assert tbl[0]["name"] == "seb skelly rocks"

    def test_string_changed_successfully_for_position(self):
        tbl = Table(names=["id", "name"], data=[[0], ["seb skelly"]])
        tbl = scopesim.utils.change_table_entry(tbl, "name", "seb skelly rocks",
                                               position=0)
        assert tbl[0]["name"] == "seb skelly rocks"

    def test_raise_error_for_neither_position_or_pattern_given(self):
        tbl = Table(names=["id", "name"], data=[[0], ["seb skelly"]])
        with pytest.raises(ValueError):
            scopesim.utils.change_table_entry(tbl, "name", "seb skelly rocks")




@pytest.mark.usefixtures("temp_directory_structure")
class TestSetUpLocalPackageDirectory:
    def test_four_folders_exist(self):
        rcnames = ["FILE_SCOPE_PKG_LOCAL_PATH", "FILE_INST_PKG_LOCAL_PATH",
                   "FILE_PSF_LOCAL_PATH",       "FILE_SRC_PKG_LOCAL_PATH"]
        for rcname in rcnames:
            filename = os.path.join(_parent_path, sim_db.rc.__rc__[rcname])
            assert os.path.exists(filename)

    def test_three_db_files_exist(self):
        for db_path in sim_db._local_paths():
            assert os.path.exists(db_path)


@pytest.mark.usefixtures("temp_directory_structure", "download_test_package")
class TestUnpackZipFile:
    def test_unzipped_folder_exists(self, download_test_package):
        local_filename = download_test_package
        sim_db.extract_package("test_package")
        assert os.path.exists(local_filename.replace(".zip", ""))

    def test_unzipped_folder_contains_files(self, download_test_package):
        local_filename = download_test_package
        sim_db.extract_package("test_package")
        num_files = len(os.listdir(local_filename.replace(".zip", "")))
        assert num_files > 0
