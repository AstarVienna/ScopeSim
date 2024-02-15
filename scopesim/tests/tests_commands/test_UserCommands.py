import os
from pathlib import Path
import pytest

from scopesim.commands.user_commands import UserCommands, patch_fake_symlinks


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestInit:
    def test_initialise_with_nothing(self):
        assert isinstance(UserCommands(), UserCommands)

    def test_initialised_when_passed_a_dict_of_properties(self):
        cmd = UserCommands(properties={"!OBS.dit": 60, "!ATMO.pwv": 9001})
        assert cmd["!ATMO.pwv"] > 9000

    def test_initialised_when_passed_a_instrument_yaml_dict(self):
        cmd = UserCommands(yamls=[{"alias": "ATMO",
                                   "properties": {"pwv": 9001}}])
        assert cmd["!ATMO.pwv"] > 9000

    def test_initialised_when_passed_a_list_of_yaml_names(self):
        cmd = UserCommands(packages=["test_package"],
                           yamls=["test_telescope.yaml"])
        assert cmd["!TEL.temperature"] > 9000

    def test_initialised_when_combining_yaml_dict_filename_properties(self):
        cmd = UserCommands(packages=["test_package"],
                           yamls=["test_telescope.yaml",
                                  {"alias": "ATMO",
                                   "properties": {"pwv": 9001}}],
                           properties={"!ATMO.pwv": 8999})
        assert cmd["!TEL.temperature"] > 9000
        assert cmd["!ATMO.pwv"] < 9000

    def test_initialise_with_correct_keywords(self):
        cmd = UserCommands(packages=["test_package"],
                           yamls=["test_instrument.yaml"],
                           properties={"!ATMO.life": 42})
        assert isinstance(cmd, UserCommands)
        assert cmd["!ATMO.life"] == 42
        assert cmd["!INST.pixel_scale"] == 0.5

    def test_initialised_with_filename_for_default_file(self):
        cmd = UserCommands(packages=["test_package"], yamls=["default.yaml"])
        assert cmd["!TEL.temperature"] < 9000
        assert len(cmd.yaml_dicts) == 7     # 3 yamls filenames + default

    def test_initialised_with_use_instrument(self):
        cmd = UserCommands(use_instrument="test_package")
        assert cmd["!TEL.temperature"] < 9000
        assert len(cmd.yaml_dicts) == 7     # 3 yamls filenames + default

    def test_mode_yamls(self):
        yamls = [{"alias": "OBS", "properties": {"modes": ["mode1"],
                                                 "life": 9001}}]
        mode_yamls = [{"name": "mode1",
                       "yamls": [{"alias": "OBS",
                                  "properties": {"life": 42}}]}]
        cmd = UserCommands(yamls=yamls, mode_yamls=mode_yamls)
        assert cmd["!OBS.life"] == 42

        print(cmd.list_modes())

    def test_throws_error_for_wrong_mode_name(self):
        with pytest.raises(ValueError):
            UserCommands(use_instrument="test_package", set_modes=["bogus"])

    def test_mode_yamls_read_from_file(self):
        cmd = UserCommands(use_instrument="test_package")
        assert cmd["!TEL.temperature"] < 9000
        assert cmd["!OBS.airmass"] == 2
        assert cmd.yaml_dicts[-1]["effects"][0]["kwargs"]["meaning_of_life"] == 42

    def test_init_through_repr(self):
        """Check whether we can recreate a UserCommand by evaluating its __repr__."""
        cmd1 = UserCommands(use_instrument="test_package")
        cmd2 = eval(repr(cmd1))
        assert cmd1 == cmd2


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestMiscFeatures:
    def test_updates_with_yaml_dict(self):
        yaml_input = {"alias": "TEL",
                      "properties": {"temperature": 8999}}
        cmd = UserCommands(use_instrument="test_package")
        cmd.update(yamls=[yaml_input])
        assert cmd["!TEL.temperature"] < 9000

    def test_update_works_via_setitem(self):
        cmd = UserCommands(use_instrument="test_package")
        cmd["!TEL.gigawatts"] = 1.21
        assert cmd["!TEL.gigawatts"] == 1.21

    def test_str(self):
        """Test whether __str__ gives a pretty result."""
        cmd = UserCommands(use_instrument="test_package")
        assert "├─" in str(cmd)


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestListLocalPackages:
    def test_all_packages_listed(self):
        from scopesim.commands import user_commands as uc2
        real_pkgs, ext_pkgs = uc2.list_local_packages(action="return")
        assert len(real_pkgs) > 0


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestTrackIpAddress:
    def test_see_if_theres_an_entry_on_the_server_log_file(self):
        _ = UserCommands(use_instrument="test_package")


def test_patch_fake_symlinks(tmp_path):
    """Setup a temporary directory with files and links."""
    # tmp_path is a fixture

    dircwd = Path.cwd()
    os.chdir(tmp_path)

    dir1 = tmp_path / "H1"
    dir1.mkdir()

    dir2 = dir1 / "H2"
    dir2.mkdir()

    # Normal file
    file1 = dir2 / "F1.txt"
    with open(file1, 'w') as f1:
        f1.write("Hello world!")

    # Empty file
    file2 = tmp_path / "F2.txt"
    with open(file2, 'w') as f2:
        f2.write("")

    # File with a line that is too long to be a link
    file3 = tmp_path / "F3.txt"
    with open(file3, 'w') as f3:
        f3.write("10 print hello; 20 goto 10" * 50)

    # A file with multiple lines
    file4 = tmp_path / "F4.txt"
    with open(file4, 'w') as f4:
        f4.write("Hello\nWorld\n")

    # With slashes. Backslashes would also work on windows,
    # but not on linux, so we just do not include that case.
    fakelink1 = tmp_path / "L1"
    with open(fakelink1, 'w') as f:
        f.write("H1/H2")

    # A real link
    reallink1 = tmp_path / "R1"
    try:
        reallink1.symlink_to(dir2)
    except OSError:
        # "A required privilege is not held by the client"
        # That is, developer mode is off.
        reallink1 = dir2

    root = list(tmp_path.parents)[-1]

    assert patch_fake_symlinks(dir1) == dir1.resolve()
    assert patch_fake_symlinks(dir2) == dir2.resolve()
    assert patch_fake_symlinks(file1) == file1.resolve()
    assert patch_fake_symlinks(file3) == file3.resolve()
    assert patch_fake_symlinks(file4) == file4.resolve()
    assert patch_fake_symlinks(fakelink1) == dir2.resolve()
    assert patch_fake_symlinks(reallink1) == dir2.resolve()
    assert patch_fake_symlinks(fakelink1 / "F1.txt") == file1.resolve()
    assert patch_fake_symlinks(reallink1 / "F1.txt") == file1.resolve()
    assert patch_fake_symlinks(root) == root.resolve()

    os.chdir(dircwd)
