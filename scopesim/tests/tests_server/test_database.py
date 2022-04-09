import os
import yaml
from tempfile import TemporaryDirectory
import numpy as np
from scopesim.server import database2 as db


def test_package_list_loads():
    pkgs = db.get_server_package_list()
    assert isinstance(pkgs, dict)
    assert "test_package" in pkgs
    assert "latest" in pkgs["test_package"]


def test_get_server_folder_contents():
    pkgs = db.get_server_folder_contents("locations")
    assert len(pkgs) > 0
    assert "Armazones" in pkgs[0]


class TestListPackages:
    def test_lists_all_packages_without_qualifier(self):
        pkgs = db.list_packages()
        assert "Armazones" in pkgs
        assert "MICADO" in pkgs

    def test_lists_only_packages_with_qualifier(self):
        pkgs = db.list_packages("Armazones")
        assert np.all(["Armazones" in pkg for pkg in pkgs])


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

    def test_downloads_github_version_of_package(self):
        pass
