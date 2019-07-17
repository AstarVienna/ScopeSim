import os
import pytest
from urllib.error import HTTPError

from scopesim import rc
from scopesim.server import database as db


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
    def test_returns_all_packages_when_nothing_specified(self):
        pkgs = db.list_packages()
        assert len(pkgs) > 0

    def test_returns_empty_list_when_url_wrong(self):
        url = rc.__config__["!SIM.file.server_base_url"][:-2]
        pkgs = db.list_packages(url)
        assert len(pkgs) == 0


class TestDownloadPackage:
    def test_downloads_package_successfully(self):
        rc.__config__["!SIM.file.local_packages_path"] = "./"
        pkg_path = "instruments/test_package.zip"

        save_path = db.download_package(pkg_path)
        assert os.path.exists(save_path)

        os.remove(save_path)
        assert not os.path.exists(save_path)

    def test_raise_error_when_package_not_found(self):
        with pytest.raises(HTTPError):
            db.download_package("instruments/bogus.zip")
