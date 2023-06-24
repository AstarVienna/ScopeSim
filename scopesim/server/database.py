# -*- coding: utf-8 -*-
"""
Functions to download instrument packages and example data
"""
import re
import logging
from datetime import date
from warnings import warn
from pathlib import Path
from typing import Optional, Union, List, Tuple, Set, Dict
# Python 3.8 doesn't yet know these things.......
# from collections.abc import Iterator, Iterable, Mapping
from typing import Iterator, Iterable, Mapping

from urllib.error import HTTPError
from urllib3.exceptions import HTTPError as HTTPError3
from more_itertools import first, last, groupby_transform

import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
import bs4

from scopesim import rc
from .github_utils import download_github_folder
from .example_data_utils import (download_example_data, list_example_data,
                                 get_server_elements)
from .download_utils import initiate_download, handle_download, handle_unzipping

_GrpVerType = Mapping[str, Iterable[str]]
_GrpItrType = Iterator[Tuple[str, List[str]]]


HTTP_RETRY_CODES = [403, 404, 429, 500, 501, 502, 503]


class ServerError(Exception):
    """Some error with the server or connection to the server."""

class PkgNotFoundError(Exception):
    """Unable to find given package or given release of that package."""

def get_server_package_list():
    warn("Function Depreciated", DeprecationWarning, stacklevel=2)

    # Emulate legacy API without using the problematic yaml file
    folders = list(dict(crawl_server_dirs()).keys())
    pkgs_dict = {}
    for dir_name in folders:
        p_list = [_parse_package_version(package) for package
                  in get_server_folder_contents(dir_name)]
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


def get_server_folder_contents(dir_name: str,
                               unique_str: str = ".zip$") -> Iterator[str]:
    url = rc.__config__["!SIM.file.server_base_url"] + dir_name

    retry_strategy = Retry(total=2,
                           status_forcelist=HTTP_RETRY_CODES,
                           allowed_methods=["GET"])
    adapter = HTTPAdapter(max_retries=retry_strategy)

    try:
        with requests.Session() as session:
            session.mount("https://", adapter)
            result = session.get(url).content
    except (requests.exceptions.ConnectionError,
            requests.exceptions.RetryError) as error:
        logging.error(error)
        raise ServerError("Cannot connect to server. "
                          f"Attempted URL was: {url}.") from error
    except Exception as error:
        logging.error(("Unhandled exception occured while accessing server."
                      "Attempted URL was: %s."), url)
        logging.error(error)
        raise error

    soup = bs4.BeautifulSoup(result, features="lxml")
    hrefs = soup.find_all("a", href=True, string=re.compile(unique_str))
    pkgs = (href.string for href in hrefs)

    return pkgs


def _get_package_name(package: str) -> str:
    return package.split(".", maxsplit=1)[0]


def _parse_raw_version(raw_version: str) -> str:
    """Catch initial package version which has no date info

    Set initial package version to basically "minus infinity".
    """
    if raw_version in ("", "zip"):
        return str(date(1, 1, 1))
    return raw_version.strip(".zip")


def _unparse_raw_version(raw_version: str, package_name: str) -> str:
    """Turn version string back into full zip folder name
    
    If initial version was set with `_parse_raw_version`, revert that.
    """
    if raw_version == str(date(1, 1, 1)):
        return f"{package_name}.zip"
    return f"{package_name}.{raw_version}.zip"


def _parse_package_version(package: str) -> Tuple[str, str]:
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


def get_all_stable(version_groups: _GrpVerType) -> Iterator[Tuple[str, str]]:
    """
    Yield the most recent version (stable or dev) of each package.

    Parameters
    ----------
    version_groups : Mapping[str, Iterable[str]]
        DESCRIPTION.

    Yields
    ------
    Iterator[Tuple[str, str]]
        Iterator of package name - latest stable version pairs.

    """
    for package_name, versions in version_groups.items():
        yield (package_name, get_stable(versions))


def get_all_latest(version_groups: _GrpVerType) -> Iterator[Tuple[str, str]]:
    """
    Yield the most recent stable (not "dev") version of each package.

    Parameters
    ----------
    version_groups : Mapping[str, Iterable[str]]
        DESCRIPTION.

    Yields
    ------
    Iterator[Tuple[str, str]]
        Iterator of package name - latest version pairs.

    """
    for package_name, versions in version_groups.items():
        yield (package_name, get_latest(versions))


def group_package_versions(all_packages: Iterable[Tuple[str, str]]) -> _GrpItrType:
    """Group different versions of packages by package name"""
    version_groups = groupby_transform(sorted(all_packages),
                                       keyfunc=first,
                                       valuefunc=last,
                                       reducefunc=list)
    return version_groups


def crawl_server_dirs() -> Iterator[Tuple[str, Set[str]]]:
    """Search all folders on server for .zip files"""
    for dir_name in get_server_folder_contents("", "/"):
        logging.info("Searching folder '%s'", dir_name)
        try:
            p_dir = get_server_folder_package_names(dir_name)
        except ValueError as err:
            logging.info(err)
            continue
        logging.info("Found packages %s.", p_dir)
        yield dir_name, p_dir


def get_all_package_versions() -> Dict[str, List[str]]:
    """Gather all versions for all packages present in any folder on server"""
    grouped = {}
    folders = list(dict(crawl_server_dirs()).keys())
    for dir_name in folders:
        p_list = [_parse_package_version(package) for package
                  in get_server_folder_contents(dir_name)]
        grouped.update(group_package_versions(p_list))
    return grouped


def get_package_folders() -> Dict[str, str]:
    folder_dict = {pkg: path.strip("/")
                   for path, pkgs in dict(crawl_server_dirs()).items()
                   for pkg in pkgs}
    return folder_dict


def get_server_folder_package_names(dir_name: str) -> Set[str]:
    """
    Retrieve all unique package names present on server in `dir_name` folder.

    Parameters
    ----------
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
                     in get_server_folder_contents(dir_name)}

    if not package_names:
        raise ValueError(f"No packages found in directory \"{dir_name}\".")

    return package_names


def get_all_packages_on_server() -> Iterator[Tuple[str, set]]:
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
    Iterator[Tuple[str, set]]
        Key-value pairs of folder and corresponding package names.

    """
    # TODO: this basically does the same as the crawl function...
    for dir_name in ("locations", "telescopes", "instruments"):
        package_names = get_server_folder_package_names(dir_name)
        yield dir_name, package_names


def list_packages(pkg_name: Optional[str] = None) -> List[str]:
    """
    List all packages, or all variants of a single package

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

    if not pkg_name in all_grouped:
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


def _download_single_package(pkg_name: str, release: str, all_versions,
                             folder_dict: Path, base_url: str, save_dir: Path,
                             padlen: int, from_cache: bool) -> Path:
    if pkg_name not in all_versions:
        raise PkgNotFoundError(f"Unable to find {release} release for "
                               f"package '{pkg_name}' on server {base_url}.")

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
    pkg_url = f"{base_url}{folder_dict[pkg_name]}/{zip_name}"

    try:
        if from_cache is None:
            from_cache = rc.__config__["!SIM.file.use_cached_downloads"]

        response = initiate_download(pkg_url, from_cache, "test_cache")
        save_path = save_dir / f"{pkg_name}.zip"
        handle_download(response, save_path, pkg_name, padlen)
        handle_unzipping(save_path, save_dir, pkg_name, padlen)

    except HTTPError3 as error:
        logging.error(error)
        msg = f"Unable to find file: {pkg_url + pkg_name}"
        raise ValueError(msg) from error
    except HTTPError as error:
        logging.error("urllib (not urllib3) error was raised, this should "
                      "not happen anymore!")
        logging.error(error)
    except requests.exceptions.ConnectionError as error:
        logging.error(error)
        raise ServerError("Cannot connect to server.") from error
    except Exception as error:
        logging.error(("Unhandled exception occured while accessing server."
                      "Attempted URL was: %s."), base_url)
        logging.error(error)
        raise error

    return save_path.absolute()


def download_packages(pkg_names: Union[Iterable[str], str],
                      release: str = "stable",
                      save_dir: Optional[str] = None,
                      from_cache: Optional[bool] = None) -> List[Path]:
    """
    Download one or more packages to the local disk

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
    base_url = rc.__config__["!SIM.file.server_base_url"]

    print("Gathering information from server ...")

    all_versions = get_all_package_versions()
    folder_dict = get_package_folders()

    print("Connection successful, starting download ...")

    if isinstance(pkg_names, str):
        pkg_names = [pkg_names]

    padlen = len(max(pkg_names, key=len))
    save_paths = []
    for pkg_name in pkg_names:
        try:
            pkg_path = _download_single_package(pkg_name, release, all_versions,
                                                folder_dict, base_url, save_dir,
                                                padlen, from_cache)
        except PkgNotFoundError as error:
            logging.error("\n")  # needed until tqdm redirect is implemented
            logging.error(error)
            logging.error("Skipping download of package '%s'", pkg_name)
            continue
        save_paths.append(pkg_path)

    return save_paths

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
    return download_packages(pkg_names, release="stable", save_dir=save_dir,
                             from_cache=from_cache)
