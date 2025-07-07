import csv
import pytest
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

from scopesim.server import ServerError
from scopesim.server import database as db
from scopesim.server import example_data_utils as dbex
from scopesim.server import github_utils as dbgh


@pytest.fixture(scope="class")
def mock_client():
    # TODO: investigate proper mocking via httpx
    with db.create_client(db.get_base_url()) as client:
        yield client


@pytest.fixture(scope="class")
def all_pkg():
    # TODO: This could use some kind of mock to avoid server access
    yield db.get_all_package_versions()


def test_package_throws():
    with pytest.raises(AttributeError):
        db.get_server_package_list()


def test_get_package_name():
    pkg_name = db._get_package_name("Packagename.2022-01-01.dev.zip")
    assert pkg_name == "Packagename"


@pytest.mark.webtest
def test_get_all_latest():
    all_pkg = db.get_all_package_versions()
    assert dict(db.get_all_latest(all_pkg))["test_package"].endswith(".dev")


@pytest.mark.webtest
class TestGetZipname:

    def test_gets_stable(self, all_pkg):
        zipname = db._get_zipname("test_package", "stable", all_pkg)
        assert zipname.startswith("test_package.")
        assert zipname.endswith(".zip")

    def test_gets_latest(self, all_pkg):
        zipname = db._get_zipname("test_package", "latest", all_pkg)
        assert zipname.startswith("test_package.")
        assert zipname.endswith(".dev.zip")

    def test_throws_for_nonexisting_release(self, all_pkg):
        with pytest.raises(ValueError):
            db._get_zipname("test_package", "bogus", all_pkg)


@pytest.mark.webtest
class TestGetServerFolderContents:
    def test_downloads_locations(self, mock_client):
        pkgs = list(db.get_server_folder_contents(mock_client, "locations"))
        assert len(pkgs) > 0

    def test_downloads_telescopes(self, mock_client):
        pkgs = list(db.get_server_folder_contents(mock_client, "telescopes"))
        assert len(pkgs) > 0

    def test_downloads_instruments(self, mock_client):
        pkgs = list(db.get_server_folder_contents(mock_client, "instruments"))
        assert len(pkgs) > 0

    def test_finds_armazones(self, mock_client):
        pkgs = list(db.get_server_folder_contents(mock_client, "locations"))
        assert "Armazones" in pkgs[0]

    def test_throws_for_wrong_url_server(self):
        with db.create_client("https://scopesim.univie.ac.at/bogus/") as client:
            with pytest.raises(ServerError):
                list(db.get_server_folder_contents(client, "locations"))


@pytest.mark.webtest
class TestListPackages:
    def test_lists_all_packages_without_qualifier(self):
        pkgs = db.list_packages()
        assert "Armazones" in pkgs
        assert "MICADO" in pkgs

    def test_lists_only_packages_with_qualifier(self):
        pkgs = db.list_packages("Armazones")
        assert all("Armazones" in pkg for pkg in pkgs)

    def test_throws_for_nonexisting_pkgname(self):
        with pytest.raises(ValueError):
            db.list_packages("bogus")


class TestDownloadPackages:
    @pytest.mark.webtest
    def test_downloads_stable_package(self):
        with TemporaryDirectory() as tmpdir:
            db.download_packages(["test_package"], release="stable",
                                 save_dir=tmpdir)
            assert Path(tmpdir, "test_package.zip").exists()

            version_path = Path(tmpdir, "test_package", "version.yaml")
            assert version_path.exists()

            with open(version_path) as f:
                version_dict = yaml.full_load(f)
            assert version_dict["release"] == "stable"

    @pytest.mark.webtest
    def test_downloads_latest_package(self):
        with TemporaryDirectory() as tmpdir:
            db.download_packages("test_package", release="latest",
                                 save_dir=tmpdir)
            version_path = Path(tmpdir, "test_package", "version.yaml")
            with version_path.open("r", encoding="utf-8") as file:
                version_dict = yaml.full_load(file)

            assert version_dict["release"] == "dev"

    @pytest.mark.webtest
    def test_downloads_specific_package(self):
        release = "2022-04-09.dev"
        with TemporaryDirectory() as tmpdir:
            db.download_packages(["test_package"], release=release,
                                 save_dir=tmpdir)
            version_path = Path(tmpdir, "test_package", "version.yaml")
            with version_path.open("r", encoding="utf-8") as file:
                version_dict = yaml.full_load(file)

            assert version_dict["version"] == release

    @pytest.mark.webtest(github=True)
    @pytest.mark.filterwarnings("ignore:Downloading IRDB packages:FutureWarning")
    @pytest.mark.filterwarnings("ignore:The function*:DeprecationWarning")
    def test_downloads_github_version_of_package_with_semicolon(self):
        release = "github:728761fc76adb548696205139e4e9a4260401dfc"
        with TemporaryDirectory() as tmpdir:
            db.download_packages("ELT", release=release,
                                 save_dir=tmpdir)
            filename = Path(tmpdir, "ELT", "EC_sky_25.tbl")

            assert filename.exists()

    @pytest.mark.webtest(github=True)
    @pytest.mark.filterwarnings("ignore:Downloading IRDB packages:FutureWarning")
    @pytest.mark.filterwarnings("ignore:The function*:DeprecationWarning")
    def test_downloads_github_version_of_package_with_at_symbol(self):
        release = "github@728761fc76adb548696205139e4e9a4260401dfc"
        with TemporaryDirectory() as tmpdir:
            db.download_packages("ELT", release=release,
                                 save_dir=tmpdir)
            filename = Path(tmpdir, "ELT", "EC_sky_25.tbl")

            assert filename.exists()


@pytest.mark.webtest(github=True)
@pytest.mark.filterwarnings("ignore:The function*:DeprecationWarning")
class TestDownloadGithubFolder:
    def test_downloads_current_package(self):
        with TemporaryDirectory() as tmpdir:
            # tmpdir = "."
            url = "https://github.com/AstarVienna/irdb/tree/dev_master/MICADO"
            dbgh.download_github_folder(url, output_dir=tmpdir)
            filename = Path(tmpdir, "MICADO", "default.yaml")

            assert filename.exists()

    def test_downloads_with_old_commit_hash(self):
        with TemporaryDirectory() as tmpdir:
            url = "https://github.com/AstarVienna/irdb/tree/728761fc76adb548696205139e4e9a4260401dfc/ELT"
            dbgh.download_github_folder(url, output_dir=tmpdir)
            filename = Path(tmpdir, "ELT", "EC_sky_25.tbl")

            assert filename.exists()

    def test_throws_for_bad_url(self):
        with TemporaryDirectory() as tmpdir:
            url = "https://github.com/AstarVienna/irdb/tree/bogus/MICADO"
            with pytest.raises(ServerError):
                dbgh.download_github_folder(url, output_dir=tmpdir)


def test_registry_files():
    registry = (Path(__file__).parent.parent.parent /
                "server/example_data_registry.txt")
    filelist = dbex.list_example_data(return_files=True)
    regfiles = dict(csv.reader(
        registry.open(encoding="utf-8"), delimiter=" ")).keys()
    assert regfiles == set(filelist)
