import pytest
import os
import sys
from tempfile import TemporaryDirectory
from urllib3.exceptions import HTTPError

import yaml
import numpy as np

from scopesim.server import database as db
from scopesim import rc


def test_package_list_loads():
    pkgs = db.get_server_package_list()
    assert isinstance(pkgs, dict)
    assert "test_package" in pkgs
    assert "latest" in pkgs["test_package"]


def test_get_server_folder_contents():
    pkgs = db.get_server_folder_contents("locations")
    assert len(pkgs) > 0
    assert "Armazones" in pkgs[0]


class TestGetServerElements:
    def test_throws_an_error_if_url_doesnt_exist(self):
        with pytest.raises(ValueError):
            db.get_server_elements(url="www.bogus.server")

    def test_returns_folders_if_server_exists(self):
        url = rc.__config__["!SIM.file.server_base_url"]
        pkgs = db.get_server_elements(url)
        assert all([loc in pkgs for loc in
                    ["locations/", "telescopes/", "instruments/"]])

    def test_returns_files_if_zips_exist(self):
        url = rc.__config__["!SIM.file.server_base_url"]
        dir = "instruments/"
        pkgs = db.get_server_elements(url + dir, ".zip")
        assert "test_package.zip" in pkgs


class TestListPackages:
    def test_lists_all_packages_without_qualifier(self):
        pkgs = db.list_packages()
        assert "Armazones" in pkgs
        assert "MICADO" in pkgs

    def test_lists_only_packages_with_qualifier(self):
        pkgs = db.list_packages("Armazones")
        assert np.all(["Armazones" in pkg for pkg in pkgs])


class TestDownloadPackage:
    """
    Old download function, for backwards compatibility
    """
    def test_downloads_package_successfully(self):
        pkg_path = "instruments/test_package.zip"

        save_paths = db.download_package(pkg_path)
        assert os.path.exists(save_paths[0])

    def test_raise_error_when_package_not_found(self):
        if sys.version_info.major >= 3:
            with pytest.raises(HTTPError):
                db.download_package("instruments/bogus.zip")


class TestDownloadPackages:
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

    def test_downloads_latest_package(self):
        with TemporaryDirectory() as tmpdir:
            db.download_packages("test_package", release="latest",
                                 save_dir=tmpdir, from_cache=False)
            version_path = os.path.join(tmpdir, "test_package", "version.yaml")
            with open(version_path) as f:
                version_dict = yaml.full_load(f)

            assert version_dict["release"] == "dev"

    def test_downloads_specific_package(self):
        release = "2022-04-09.dev"
        with TemporaryDirectory() as tmpdir:
            db.download_packages(["test_package"], release=release,
                                 save_dir=tmpdir, from_cache=False)
            version_path = os.path.join(tmpdir, "test_package", "version.yaml")
            with open(version_path) as f:
                version_dict = yaml.full_load(f)

            assert version_dict["version"] == release

    def test_downloads_github_version_of_package_with_semicolon(self):
        release = "github:728761fc76adb548696205139e4e9a4260401dfc"
        with TemporaryDirectory() as tmpdir:
            db.download_packages("ELT", release=release,
                                 save_dir=tmpdir, from_cache=False)
            filename = os.path.join(tmpdir, "ELT", "EC_sky_25.tbl")

            assert os.path.exists(filename)

    def test_downloads_github_version_of_package_with_at_symbol(self):
        release = "github@728761fc76adb548696205139e4e9a4260401dfc"
        with TemporaryDirectory() as tmpdir:
            db.download_packages("ELT", release=release,
                                 save_dir=tmpdir, from_cache=False)
            filename = os.path.join(tmpdir, "ELT", "EC_sky_25.tbl")

            assert os.path.exists(filename)


class TestDownloadGithubFolder:
    def test_downloads_current_package(self):
        with TemporaryDirectory() as tmpdir:
            # tmpdir = "."
            url = "https://github.com/AstarVienna/irdb/tree/dev_master/MICADO"
            db.download_github_folder(url, output_dir=tmpdir)
            filename = os.path.join(tmpdir, "MICADO", "default.yaml")

            assert os.path.exists(filename)

    def test_downloads_with_old_commit_hash(self):
        with TemporaryDirectory() as tmpdir:
            url = "https://github.com/AstarVienna/irdb/tree/728761fc76adb548696205139e4e9a4260401dfc/ELT"
            db.download_github_folder(url, output_dir=tmpdir)
            filename = os.path.join(tmpdir, "ELT", "EC_sky_25.tbl")

            assert os.path.exists(filename)


def test_old_download_package_signature():
    with TemporaryDirectory() as tmpdir:
        db.download_package(["instruments/test_package.zip"], save_dir=tmpdir)
        version_path = os.path.join(tmpdir, "test_package", "version.yaml")
        with open(version_path) as f:
            version_dict = yaml.full_load(f)

        assert version_dict["release"] == "stable"
