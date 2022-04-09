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


def test_list_packages():
    print(db.list_packages())