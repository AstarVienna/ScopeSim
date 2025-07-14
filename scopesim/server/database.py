# -*- coding: utf-8 -*-
"""Functions to download instrument packages and example data."""

from datetime import date
from warnings import warn
from pathlib import Path
from typing import Optional, Union
from collections.abc import Iterator, Iterable, Mapping

import yaml
from more_itertools import first, last, groupby_transform

from scopesim import rc
from .github_utils import download_github_folder
from .example_data_utils import download_example_data, list_example_data
from .download_utils import (get_server_folder_contents, handle_download,
                             handle_unzipping, create_client, ServerError)
from ..utils import get_logger
from ..commands.user_commands import patch_fake_symlinks


logger = get_logger(__name__)

_GrpVerType = Mapping[str, Iterable[str]]
_GrpItrType = Iterator[tuple[str, list[str]]]


class PkgNotFoundError(Exception):
    """Unable to find given package or given release of that package."""


def get_base_url():
    """Get instrument package server URL from rc.__config__."""
    return rc.__config__["!SIM.file.server_base_url"]


def get_server_package_list():
    warn("Function Depreciated", DeprecationWarning, stacklevel=2)

    # Emulate legacy API without using the problematic yaml file
    with create_client(get_base_url()) as client:
        folders = list(dict(crawl_server_dirs(client)).keys())
        pkgs_dict = {}
        for dir_name in folders:
            p_list = [_parse_package_version(package) for package
                      in get_server_folder_contents(client, dir_name)]
            grouped = dict(group_package_versions(p_list))
            for p_name in grouped:
                p_dict = {
                    "latest": _unparse_raw_version(get_latest(grouped[p_name]),
                                                   p_name).strip(".zip"),
                    "path": dir_name.strip("/"),
                    "stable": _unparse_raw_version(get_stable(grouped[p_name]),
                                                   p_name).strip(".zip"),
                }
                pkgs_dict[p_name] = p_dict

    return pkgs_dict


def _get_package_name(package: str) -> str:
    return package.split(".", maxsplit=1)[0]


def _parse_raw_version(raw_version: str) -> str:
    """Catch initial package version which has no date info.

    Set initial package version to basically "minus infinity".
    """
    if raw_version in ("", "zip"):
        return str(date(1, 1, 1))
    return raw_version.strip(".zip")


def _unparse_raw_version(raw_version: str, package_name: str) -> str:
    """Turn version string back into full zip folder name.

    If initial version was set with `_parse_raw_version`, revert that.
    """
    if raw_version == str(date(1, 1, 1)):
        return f"{package_name}.zip"
    return f"{package_name}.{raw_version}.zip"


def _parse_package_version(package: str) -> tuple[str, str]:
    p_name, p_version = package.split(".", maxsplit=1)
    return p_name, _parse_raw_version(p_version)


def _is_stable(package_version: str) -> bool:
    return not package_version.endswith("dev")


def get_stable(versions: Iterable[str]) -> str:
    """Return the most recent stable (not "dev") version."""
    return max(version for version in versions if _is_stable(version))


def get_latest(versions: Iterable[str]) -> str:
    """Return the most recent version (stable or dev)."""
    return max(versions)


def get_all_stable(version_groups: _GrpVerType) -> Iterator[tuple[str, str]]:
    """
    Yield the most recent version (stable or dev) of each package.

    Parameters
    ----------
    version_groups : Mapping[str, Iterable[str]]
        DESCRIPTION.

    Yields
    ------
    Iterator[tuple[str, str]]
        Iterator of package name - latest stable version pairs.

    """
    for package_name, versions in version_groups.items():
        yield (package_name, get_stable(versions))


def get_all_latest(version_groups: _GrpVerType) -> Iterator[tuple[str, str]]:
    """
    Yield the most recent stable (not "dev") version of each package.

    Parameters
    ----------
    version_groups : Mapping[str, Iterable[str]]
        DESCRIPTION.

    Yields
    ------
    Iterator[tuple[str, str]]
        Iterator of package name - latest version pairs.

    """
    for package_name, versions in version_groups.items():
        yield (package_name, get_latest(versions))


def group_package_versions(all_packages: Iterable[tuple[str, str]]) -> _GrpItrType:
    """Group different versions of packages by package name."""
    version_groups = groupby_transform(sorted(all_packages),
                                       keyfunc=first,
                                       valuefunc=last,
                                       reducefunc=list)
    return version_groups


def crawl_server_dirs(client=None) -> Iterator[tuple[str, set[str]]]:
    """Search all folders on server for .zip files."""
    if client is None:
        with create_client(get_base_url()) as client:
            yield from crawl_server_dirs(client)
        return

    for dir_name in get_server_folder_contents(client, "", "/"):
        logger.debug("Searching folder '%s'", dir_name)
        try:
            p_dir = get_server_folder_package_names(client, dir_name)
        except ValueError as err:
            logger.debug(err)
            continue
        logger.debug("Found packages %s.", p_dir)
        yield dir_name, p_dir


def get_all_package_versions(client=None) -> dict[str, list[str]]:
    """Gather all versions for all packages present in any folder on server."""
    if client is None:
        with create_client(get_base_url()) as client:
            return get_all_package_versions(client)

    grouped = {}
    folders = list(dict(crawl_server_dirs(client)).keys())
    for dir_name in folders:
        p_list = [_parse_package_version(package) for package
                  in get_server_folder_contents(client, dir_name)]
        grouped.update(group_package_versions(p_list))
    return grouped


def get_package_folders(client) -> dict[str, str]:
    """Map package names to server locations."""
    folders_dict = {pkg: path.strip("/")
                    for path, pkgs in dict(crawl_server_dirs(client)).items()
                    for pkg in pkgs}
    return folders_dict


def get_server_folder_package_names(client, dir_name: str) -> set[str]:
    """
    Retrieve all unique package names present on server in `dir_name` folder.

    Parameters
    ----------
    client : httpx.Client
        Pre-existing httpx Client context manager.

    dir_name : str
        Name of the folder on the server.

    Raises
    ------
    ValueError
        Raised if no valid packages are found in the given folder.

    Returns
    -------
    package_names : set of str
        Set of unique package names in `dir_name` folder.

    """
    package_names = {package.split(".", maxsplit=1)[0] for package
                     in get_server_folder_contents(client, dir_name)}

    if not package_names:
        raise ValueError(f"No packages found in directory \"{dir_name}\".")

    return package_names


def get_all_packages_on_server() -> Iterator[tuple[str, set]]:
    """
    Retrieve all unique package names present on server in known folders.

    Currently hardcoded to look in folders "locations", "telescopes" and
    "instruments". Any packages not in these folders are not returned.

    This generator function yields key-value pairs, containing the folder name
    as the key and the set of unique package names in value. Recommended useage
    is to turn the generator into a dictionary, i.e.:

    ::
        package_dict = dict(get_all_packages_on_server())

    Yields
    ------
    Iterator[tuple[str, set]]
        Key-value pairs of folder and corresponding package names.

    """
    yield from crawl_server_dirs()


def list_packages(pkg_name: Optional[str] = None) -> list[str]:
    """
    List all packages, or all variants of a single package.

    Parameters
    ----------
    pkg_name : str, optional
        - None: lists all stable packages on the server
        - <PackageName>: lists all variants of <PackageName> on the server

    Returns
    -------
    pkg_names : list

    Examples
    --------
    ::
        from scopesim import list_packages

        # list all stable packages on the server
        list_packages()

        # list all variants of a specific package
        list_packages("Armazones")

    """
    all_grouped = get_all_package_versions()

    if pkg_name is None:
        # Return all packages with any stable version
        all_stable = list(dict(get_all_stable(all_grouped)).keys())
        return all_stable

    if pkg_name not in all_grouped:
        raise ValueError(f"Package name {pkg_name} not found on server.")

    p_versions = [_unparse_raw_version(version, pkg_name)
                  for version in all_grouped[pkg_name]]
    return p_versions


def _get_zipname(pkg_name: str, release: str, all_versions) -> str:
    if release == "stable":
        zip_name = get_stable(all_versions[pkg_name])
    elif release == "latest":
        zip_name = get_latest(all_versions[pkg_name])
    else:
        release = _parse_raw_version(release)
        if release not in all_versions[pkg_name]:
            msg = (f"Requested version '{release}' of '{pkg_name}' package"
                   " could not be found on the server. Available versions "
                   f"are: {all_versions[pkg_name]}")
            raise ValueError(msg)
        zip_name = release
    return _unparse_raw_version(zip_name, pkg_name)


def _download_single_package(client, pkg_name: str, release: str, all_versions,
                             folders_dict: Path, save_dir: Path,
                             padlen: int) -> Path:
    if pkg_name not in all_versions:
        maybe = ""
        for key in folders_dict:
            if pkg_name in key or key in pkg_name:
                maybe = f"\nDid you mean '{key}' instead of '{pkg_name}'?"

        raise PkgNotFoundError(
            f"Unable to find {release} release for package '{pkg_name}' on "
            f"server {client.base_url!s}.{maybe}")

    if save_dir is None:
        save_dir = rc.__config__["!SIM.file.local_packages_path"]
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    if "github" in release:
        base_url = "https://github.com/AstarVienna/irdb/tree/"
        github_hash = release.split(":")[-1].split("@")[-1]
        pkg_url = f"{base_url}{github_hash}/{pkg_name}"
        download_github_folder(repo_url=pkg_url, output_dir=save_dir)
        return save_dir.absolute()

    zip_name = _get_zipname(pkg_name, release, all_versions)
    pkg_url = f"{folders_dict[pkg_name]}/{zip_name}"

    save_path = save_dir / f"{pkg_name}.zip"
    handle_download(client, pkg_url, save_path, pkg_name, padlen)
    handle_unzipping(save_path, save_dir, pkg_name, padlen)

    return save_path.absolute()


def download_packages(pkg_names: Union[Iterable[str], str],
                      release: str = "stable",
                      save_dir: Optional[str] = None) -> list[Path]:
    """
    Download one or more packages to the local disk.

    1. Download stable, dev
    2. Download specific version
    3. Download from github via url

    Parameters
    ----------
    pkg_names : str, list
        A list of package names, see ``list_packages()``

    release : str, optional
        By default, the most recent stable version of a package is downloaded.
        Other options are:
        - "stable" : the most recent stable version
        - "latest" : the latest development version to be published
        - a specific package filename as given by list_packages (see examples)
        - a github url for the specific branch and package (see examples)

    save_dir : str, optional
        The place on the local disk where the ``.zip`` package is to be saved.
        If left as None, defaults to the value in
        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    from_cache : bool, optional
        Use the cached versions of the packages. If None, defaults to the RC
        value: ``!SIM.file.use_cached_downloads``

    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

    Examples
    --------
    ::
        from scopesim import download_packages, list_packages

        # Stable release of a list of packages
        download_packages(["test_package", "test_package"])

        # Development release of a single package
        download_packages("test_package", release="latest")

        # Specific version of the package
        list_packages("test_package")
        download_packages("test_package", release="2022-04-09.dev")

        # Specific package from a Gtihub commit hash or branch/tag name (use "@" or ":")
        download_packages("ELT", release="github:728761fc76adb548696205139e4e9a4260401dfc")
        download_packages("ELT", release="github@728761fc76adb548696205139e4e9a4260401dfc")
        download_packages("ELT", release="github@dev_master")

    """
    base_url = get_base_url()
    logger.info("Gathering information from server ...")
    logger.debug("Accessing %s", base_url)

    with create_client(base_url) as client:
        all_versions = get_all_package_versions(client)
        folders_dict = get_package_folders(client)

        logger.info("Connection successful, starting download ...")

        if isinstance(pkg_names, str):
            pkg_names = [pkg_names]

        padlen = max(len(name) for name in pkg_names)
        save_paths = []
        for pkg_name in pkg_names:
            try:
                pkg_path = _download_single_package(
                    client, pkg_name, release, all_versions, folders_dict,
                    save_dir, padlen)
            except PkgNotFoundError as error:
                # Whole stack trace not useful for enduser.
                # Could log it to file though...
                # logger.exception(error)
                logger.error(error)
                logger.warning("Skipping download of package '%s'", pkg_name)
                continue
            save_paths.append(pkg_path)

    return save_paths


# TODO: Look into just using this in download_packages as well...
def download_missing_pkgs(instrument: str) -> None:
    """Download instrument package and required support packages."""
    # First download package itself
    zip_file = download_packages([instrument])[0]

    # If package needs other packages, download them as well
    defyam = zip_file.with_suffix("") / "default.yaml"
    with defyam.open() as file:
        pkgs = next(yaml.load_all(file, yaml.SafeLoader))["packages"]
    pkgs.remove(instrument)

    # Only actually download if necessary
    if pkgs:
        download_packages(pkgs)


def check_packages(instrument: str, download_missing: bool) -> None:
    """Check if required package is in CWD, download if needed."""
    pkgdir = Path(rc.__config__["!SIM.file.local_packages_path"])
    if not pkgdir.exists():
        pkgdir.mkdir()
    pkgdir = patch_fake_symlinks(pkgdir)
    if not (pkgdir / instrument).exists():
        logger.warning("IRDB package for %s not found.", instrument)
        if download_missing:
            download_missing_pkgs(instrument)
        else:
            raise ValueError(
                f"IRDB package for {instrument} not found, auto-download is "
                "disabled. Please set package directory or download package.")


# ==============================================================================
# Funtions below from from OLD_database.py
# ==============================================================================

# for backwards compatibility
def download_package(pkg_path, save_dir=None, url=None, from_cache=None):
    """
    DEPRECATED -- only kept for backwards compatibility

    Downloads a package to the local disk

    Parameters
    ----------
    pkg_path : str, list
        A ``.zip`` package path as given by ``list_packages()``

    save_dir : str
        The place on the local disk where the ``.zip`` package is to be saved.
        If left as None, defaults to the value in
        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    url : str
        The URL of the IRDB HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    from_cache : bool
        Use the cached versions of the packages. If None, defaults to the RC
        value: ``!SIM.file.use_cached_downloads``

    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

    """
    warn("Function Depreciated --> please use scopesim.download_package-s-()",
         DeprecationWarning, stacklevel=2)

    if isinstance(pkg_path, str):
        pkg_path = [pkg_path]

    pkg_names = [pkg.replace(".zip", "").split("/")[-1] for pkg in pkg_path]
    return download_packages(pkg_names, release="stable", save_dir=save_dir)
