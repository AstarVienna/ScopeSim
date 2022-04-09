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

