import pytest

from scopesim import rc
from scopesim.server import database2 as db


class TestGetServerElements:
    def test_returns_nothing_if_server_doesnt_respond(self):
        pkgs = db.get_server_elements(url=None)
        assert len(pkgs) == 0 and isinstance(pkgs, list)

    def test_returns_folders_if_server_exists(self):
        url = rc.__rc__["FILE_SERVER_BASE_URL"]
        pkgs = db.get_server_elements(url)
        assert "telescopes/" in pkgs
        assert "instruments/" in pkgs

    def test_returns_files_if_zips_exist(self):
        url = rc.__rc__["FILE_SERVER_BASE_URL"]
        dir = rc.__rc__["FILE_INST_PKG_LOCAL_PATH"]
        pkgs = db.get_server_elements(url + dir, ".zip")
        assert "test_package.zip" in pkgs


class TestGetPackages:
    def test_returns_all_packages_when_nothing_specified(self):
        print(db.get_packages())


