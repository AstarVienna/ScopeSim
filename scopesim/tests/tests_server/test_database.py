import pytest
import os
from tempfile import TemporaryDirectory
from urllib3.exceptions import HTTPError

import yaml
import numpy as np

from scopesim.server import database as db
from scopesim.server import example_data_utils as dbex
from scopesim.server import github_utils as dbgh
from scopesim import rc


@pytest.mark.webtest
def test_package_list_loads():
    pkgs = db.get_server_package_list()
    assert isinstance(pkgs, dict)
    assert "test_package" in pkgs
    assert "latest" in pkgs["test_package"]


def test_get_package_name():
    pkg_name = db._get_package_name("Packagename.2022-01-01.dev.zip")
    assert pkg_name == "Packagename"


@pytest.mark.webtest
def test_get_all_latest():
    all_pkg = db.get_all_package_versions()
    assert dict(db.get_all_latest(all_pkg))["test_package"].endswith(".dev")


@pytest.mark.webtest
class TestGetZipname:
    # TODO: This could use some kind of mock to avoid server access
    all_pkg = db.get_all_package_versions()

    def test_gets_stable(self):
        zipname = db._get_zipname("test_package", "stable", self.all_pkg)
        assert zipname.startswith("test_package.")
        assert zipname.endswith(".zip")

    def test_gets_latest(self):
        zipname = db._get_zipname("test_package", "latest", self.all_pkg)
        assert zipname.startswith("test_package.")
        assert zipname.endswith(".dev.zip")

    def test_throws_for_nonexisting_release(self):
        with pytest.raises(ValueError):
            db._get_zipname("test_package", "bogus", self.all_pkg)


class TestGetServerFolderContents:
    @pytest.mark.webtest
    def test_downloads_locations(self):
        pkgs = list(db.get_server_folder_contents("locations"))
        assert len(pkgs) > 0

    @pytest.mark.webtest
    def test_downloads_telescopes(self):
        pkgs = list(db.get_server_folder_contents("telescopes"))
        assert len(pkgs) > 0

    @pytest.mark.webtest
    def test_downloads_instruments(self):
        pkgs = list(db.get_server_folder_contents("instruments"))
        assert len(pkgs) > 0

    @pytest.mark.webtest
    def test_finds_armazones(self):
        pkgs = list(db.get_server_folder_contents("locations"))
        assert "Armazones" in pkgs[0]

    @pytest.mark.webtest
    def test_throws_for_wrong_url_server(self):
        original_url = rc.__config__["!SIM.file.server_base_url"]
        rc.__config__["!SIM.file.server_base_url"] = "https://scopesim.univie.ac.at/bogus/"
        with pytest.raises(db.ServerError):
            list(db.get_server_folder_contents("locations"))
        rc.__config__["!SIM.file.server_base_url"] = original_url


class TestGetServerElements:
    @pytest.mark.webtest
    def test_throws_an_error_if_url_doesnt_exist(self):
        with pytest.raises(ValueError):
            dbex.get_server_elements(url="www.bogus.server")

    @pytest.mark.webtest
    def test_returns_folders_if_server_exists(self):
        url = rc.__config__["!SIM.file.server_base_url"]
        pkgs = dbex.get_server_elements(url)
        assert all([loc in pkgs for loc in
                    ["locations/", "telescopes/", "instruments/"]])

    @pytest.mark.webtest
    def test_returns_files_if_zips_exist(self):
        url = rc.__config__["!SIM.file.server_base_url"]
        dir = "instruments/"
        pkgs = dbex.get_server_elements(url + dir, ".zip")
        assert "test_package.zip" in pkgs


class TestListPackages:
    @pytest.mark.webtest
    def test_lists_all_packages_without_qualifier(self):
        pkgs = db.list_packages()
        assert "Armazones" in pkgs
        assert "MICADO" in pkgs

    @pytest.mark.webtest
    def test_lists_only_packages_with_qualifier(self):
        pkgs = db.list_packages("Armazones")
        assert np.all(["Armazones" in pkg for pkg in pkgs])

    @pytest.mark.webtest
    def test_throws_for_nonexisting_pkgname(self):
        with pytest.raises(ValueError):
            db.list_packages("bogus")


class TestDownloadPackage:
    """
    Old download function, for backwards compatibility
    """
    @pytest.mark.webtest
    def test_downloads_package_successfully(self):
        pkg_path = "instruments/test_package.zip"

        save_paths = db.download_package(pkg_path)
        assert os.path.exists(save_paths[0])

    # This no longer raises, but logs an error. This is intended.
    # TODO: Change test to capture log and assert if error log is present.
    # Actually, the new single download function should be tested here instead
    # def test_raise_error_when_package_not_found(self):
    #     if sys.version_info.major >= 3:
    #         with pytest.raises(HTTPError):
    #             db.download_package("instruments/bogus.zip")


class TestDownloadPackages:
    @pytest.mark.webtest
    def test_downloads_stable_package(self):
        with TemporaryDirectory() as tmpdir:
            db.download_packages(["test_package"], release="stable",
                                 save_dir=tmpdir, from_cache=False)
            assert os.path.exists(os.path.join(tmpdir, "test_package.zip"))

            version_path = os.path.join(tmpdir, "test_package", "version.yaml")
            assert os.path.exists(version_path)

            with open(version_path) as f:
                version_dict = yaml.full_load(f)
            assert version_dict["release"] == "stable"

    @pytest.mark.webtest
    def test_downloads_latest_package(self):
        with TemporaryDirectory() as tmpdir:
            db.download_packages("test_package", release="latest",
                                 save_dir=tmpdir, from_cache=False)
            version_path = os.path.join(tmpdir, "test_package", "version.yaml")
            with open(version_path) as f:
                version_dict = yaml.full_load(f)

            assert version_dict["release"] == "dev"

    @pytest.mark.webtest
    def test_downloads_specific_package(self):
        release = "2022-04-09.dev"
        with TemporaryDirectory() as tmpdir:
            db.download_packages(["test_package"], release=release,
                                 save_dir=tmpdir, from_cache=False)
            version_path = os.path.join(tmpdir, "test_package", "version.yaml")
            with open(version_path) as f:
                version_dict = yaml.full_load(f)

            assert version_dict["version"] == release

    @pytest.mark.skip(reason="fails too often with timeout")
    @pytest.mark.webtest
    def test_downloads_github_version_of_package_with_semicolon(self):
        release = "github:728761fc76adb548696205139e4e9a4260401dfc"
        with TemporaryDirectory() as tmpdir:
            db.download_packages("ELT", release=release,
                                 save_dir=tmpdir, from_cache=False)
            filename = os.path.join(tmpdir, "ELT", "EC_sky_25.tbl")

            assert os.path.exists(filename)

    @pytest.mark.skip(reason="fails too often with timeout")
    @pytest.mark.webtest
    def test_downloads_github_version_of_package_with_at_symbol(self):
        release = "github@728761fc76adb548696205139e4e9a4260401dfc"
        with TemporaryDirectory() as tmpdir:
            db.download_packages("ELT", release=release,
                                 save_dir=tmpdir, from_cache=False)
            filename = os.path.join(tmpdir, "ELT", "EC_sky_25.tbl")

            assert os.path.exists(filename)


@pytest.mark.skip(reason="fails too often with timeout")
class TestDownloadGithubFolder:
    @pytest.mark.webtest
    def test_downloads_current_package(self):
        with TemporaryDirectory() as tmpdir:
            # tmpdir = "."
            url = "https://github.com/AstarVienna/irdb/tree/dev_master/MICADO"
            dbgh.download_github_folder(url, output_dir=tmpdir)
            filename = os.path.join(tmpdir, "MICADO", "default.yaml")

            assert os.path.exists(filename)

    @pytest.mark.webtest
    def test_downloads_with_old_commit_hash(self):
        with TemporaryDirectory() as tmpdir:
            url = "https://github.com/AstarVienna/irdb/tree/728761fc76adb548696205139e4e9a4260401dfc/ELT"
            dbgh.download_github_folder(url, output_dir=tmpdir)
            filename = os.path.join(tmpdir, "ELT", "EC_sky_25.tbl")

            assert os.path.exists(filename)

    @pytest.mark.webtest
    def test_throws_for_bad_url(self):
        with TemporaryDirectory() as tmpdir:
            url = "https://github.com/AstarVienna/irdb/tree/bogus/MICADO"
            with pytest.raises(dbgh.ServerError):
                dbgh.download_github_folder(url, output_dir=tmpdir)


@pytest.mark.webtest
def test_old_download_package_signature():
    with TemporaryDirectory() as tmpdir:
        db.download_package(["instruments/test_package.zip"], save_dir=tmpdir)
        version_path = os.path.join(tmpdir, "test_package", "version.yaml")
        with open(version_path) as f:
            version_dict = yaml.full_load(f)

        assert version_dict["release"] == "stable"
