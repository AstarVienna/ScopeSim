"""
Functions to download instrument packages and example data
"""
import json
import re
import shutil
import urllib.request
from zipfile import ZipFile
import logging
from datetime import date
from warnings import warn
from pathlib import Path
from typing import Optional, Union, List, Tuple, Set, Dict
from collections.abc import Iterator, Iterable, Mapping
from shutil import get_terminal_size

from urllib.error import HTTPError
from urllib3.exceptions import HTTPError as HTTPError3
from more_itertools import first, last, groupby_transform

import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from requests_cache import CachedSession
import bs4
from astropy.utils.data import download_file
from tqdm import tqdm
# from tqdm.contrib.logging import logging_redirect_tqdm

from scopesim import rc


GrpVerType = Mapping[str, Iterable[str]]
GrpItrType = Iterator[Tuple[str, List[str]]]


width, _ = get_terminal_size((50, 20))
width = int(.8 * width)
bar_width = max(width - 50, 10)
tqdm_fmt = f"{{l_bar}}{{bar:{bar_width}}}{{r_bar}}{{bar:-{bar_width}b}}"


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

    try:
        result = requests.get(url).content
    except Exception as error:
        raise ValueError(f"URL returned error: {url}") from error

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


def get_all_stable(version_groups: GrpVerType) -> Iterator[Tuple[str, str]]:
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


def get_all_latest(version_groups: GrpVerType) -> Iterator[Tuple[str, str]]:
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


def group_package_versions(all_packages: Iterable[Tuple[str, str]]) -> GrpItrType:
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


def _create_session(cached: bool = False, cache_name: str = ""):
    if cached:
        return CachedSession(cache_name)
    return requests.Session()


def _handle_download(response, save_path: Path, pkg_name: str,
                     padlen: int, chunk_size: int = 128) -> None:
    try:
        with tqdm.wrapattr(save_path.open("wb"), "write", miniters=1,
                           total=int(response.headers.get("content-length", 0)),
                           desc=f"Downloading {pkg_name:<{padlen}}",
                           bar_format=tqdm_fmt) as file:
            for chunk in response.iter_content(chunk_size=chunk_size):
                file.write(chunk)
    except Exception as error:
        raise error
    finally:
        file.close()


def _handle_unzipping(save_path: Path, pkg_name: str, padlen: int) -> None:
    with ZipFile(save_path, "r") as zip_ref:
        namelist = zip_ref.namelist()
        for file in tqdm(iterable=namelist, total=len(namelist),
                         desc=f"Extracting  {pkg_name:<{padlen}}",
                         bar_format=tqdm_fmt):
            zip_ref.extract(member=file)


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

    all_versions = get_all_package_versions()
    folder_dict = get_package_folders()

    if isinstance(pkg_names, str):
        pkg_names = [pkg_names]

    padlen = len(max(pkg_names, key=len))
    save_paths = []
    for pkg_name in pkg_names:
        if pkg_name not in all_versions:
            raise HTTPError3(f"Unable to find {release} release for package: "
                             f"{base_url + pkg_name}")

        if save_dir is None:
            save_dir = rc.__config__["!SIM.file.local_packages_path"]
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

        if "github" in release:
            base_url = "https://github.com/AstarVienna/irdb/tree/"
            github_hash = release.split(":")[-1].split("@")[-1]
            pkg_url = f"{base_url}{github_hash}/{pkg_name}"
            download_github_folder(repo_url=pkg_url, output_dir=save_dir)
            save_paths.append(save_dir.absolute())
            continue

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

        zip_name = _unparse_raw_version(zip_name, pkg_name)
        pkg_url = f"{base_url}{folder_dict[pkg_name]}/{zip_name}"

        try:
            if from_cache is None:
                from_cache = rc.__config__["!SIM.file.use_cached_downloads"]

            retry_strategy = Retry(total=5, backoff_factor=2,
                                   status_forcelist=[429, 500, 501, 502, 503],
                                   allowed_methods=["GET"])
            adapter = HTTPAdapter(max_retries=retry_strategy)
            with _create_session(from_cache, "foo") as session:
                session.mount("https://", adapter)
                response = session.get(pkg_url, stream=True)

            save_path = save_dir / f"{pkg_name}.zip"
            _handle_download(response, save_path, pkg_name, padlen)
            _handle_unzipping(save_path, pkg_name, padlen)

        except (HTTPError, HTTPError3) as error:
            msg = f"Unable to find file: {pkg_url + pkg_name}"
            raise ValueError(msg) from error

        save_paths.append(save_path.absolute())

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

def get_server_elements(url: str, unique_str: str = "/") -> List[str]:
    """
    Returns a list of file and/or directory paths on the HTTP server ``url``

    Parameters
    ----------
    url : str
        The URL of the IRDB HTTP server.

    unique_str : str, list
        A unique string to look for in the beautiful HTML soup:
        "/" for directories this, ".zip" for packages

    Returns
    -------
    paths : list
        List of paths containing in ``url`` which contain ``unique_str``

    """
    if isinstance(unique_str, str):
        unique_str = [unique_str]

    try:
        result = requests.get(url).content
    except Exception as error:
        raise ValueError(f"URL returned error: {url}") from error

    soup = bs4.BeautifulSoup(result, features="lxml")
    paths = soup.findAll("a", href=True)
    select_paths = []
    for the_str in unique_str:
        select_paths += [tmp.string for tmp in paths
                         if tmp.string is not None and the_str in tmp.string]
    return select_paths


def list_example_data(url: Optional[str] = None,
                      return_files: bool = False,
                      silent: bool = False) -> List[str]:
    """
    List all example files found under ``url``

    Parameters
    ----------
    url : str
        The URL of the database HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    return_files : bool
        If True, returns a list of file names

    silent : bool
        If True, does not print the list of file names

    Returns
    -------
    all_files : list of str
        A list of paths to the example files relative to ``url``.
        The full string should be passed to ``download_example_data``.
    """

    def print_file_list(the_files, loc=""):
        print(f"\nFiles saved {loc}\n" + "=" * (len(loc) + 12))
        for _file in the_files:
            print(_file)

    if url is None:
        url = rc.__config__["!SIM.file.server_base_url"]

    return_file_list = []
    server_files = []
    folders = get_server_elements(url, "example_data")
    for folder in folders:
        files = get_server_elements(url + folder, ("fits", "txt", "dat"))
        server_files += files
    if not silent:
        print_file_list(server_files, f"on the server: {url + 'example_data/'}")
        return_file_list += server_files

    if return_files:
        return return_file_list

    return None


def download_example_data(file_path: Union[Iterable[str], str],
                          save_dir: Optional[Union[Path, str]] = None,
                          url: Optional[str] = None,
                          from_cache: Optional[bool] = None) -> List[Path]:
    """
    Downloads example fits files to the local disk

    Parameters
    ----------
    file_path : str, list
        Name(s) of FITS file(s) as given by ``list_example_data()``

    save_dir : str
        The place on the local disk where the downloaded files are to be saved.
        If left as None, defaults to the current working directory.

    url : str
        The URL of the database HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    from_cache : bool
        Use the cached versions of the files. If None, defaults to the RC
        value: ``!SIM.file.use_cached_downloads``

    Returns
    -------
    save_path : Path or list of Paths
        The absolute path(s) to the saved files
    """
    if isinstance(file_path, Iterable) and not isinstance(file_path, str):
        # Recursive
        save_path = [download_example_data(thefile, save_dir, url)
                     for thefile in file_path]
        return save_path

    if not isinstance(file_path, str):
        raise TypeError("file_path must be str or iterable of str, found "
                        f"{type(file_path) = }")

    if url is None:
        url = rc.__config__["!SIM.file.server_base_url"]
    if save_dir is None:
        save_dir = Path.cwd()
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    file_path = Path(file_path)

    try:
        if from_cache is None:
            from_cache = rc.__config__["!SIM.file.use_cached_downloads"]
        cache_path = download_file(f"{url}example_data/{file_path}",
                                   cache=from_cache)
        save_path = save_dir / file_path.name
        file_path = shutil.copy2(cache_path, str(save_path))
    except (HTTPError, HTTPError3) as error:
        msg = f"Unable to find file: {url + 'example_data/' + file_path}"
        raise ValueError(msg) from error

    save_path = save_path.absolute()
    return save_path


# """
# 2022-04-10 (KL)
# Code taken directly from https://github.com/sdushantha/gitdir
# Adapted for ScopeSim usage.
# Many thanks to the authors!
# """

def create_github_url(url: str) -> None:
    """
    From the given url, produce a URL that is compatible with Github's REST API.

    Can handle blob or tree paths.
    """
    repo_only_url = re.compile(r"https:\/\/github\.com\/[a-z\d](?:[a-z\d]|-(?=[a-z\d])){0,38}\/[a-zA-Z0-9]+$")
    re_branch = re.compile("/(tree|blob)/(.+?)/")

    # Check if the given url is a url to a GitHub repo. If it is, tell the
    # user to use 'git clone' to download it
    if re.match(repo_only_url,url):
        message = ("✘ The given url is a complete repository. Use 'git clone'"
                   " to download the repository")
        logging.error(message)
        raise ValueError(message)

    # extract the branch name from the given url (e.g master)
    branch = re_branch.search(url)
    download_dirs = url[branch.end():]
    api_url = (url[:branch.start()].replace("github.com", "api.github.com/repos", 1) +
               f"/contents/{download_dirs}?ref={branch.group(2)}")
    return api_url, download_dirs


def download_github_folder(repo_url: str,
                           output_dir: Union[Path, str] = "./") -> None:
    """
    Downloads the files and directories in repo_url.

    Re-written based on the on the download function
    `here <https://github.com/sdushantha/gitdir/blob/f47ce9d85ee29f8612ce5ae804560a12b803ddf3/gitdir/gitdir.py#L55>`_
    """
    output_dir = Path(output_dir)

    # convert repo_url into an api_url
    api_url, download_dirs = create_github_url(repo_url)

    # get the contents of the github folder
    user_interrupt_text = "GitHub download interrupted by User"
    try:
        opener = urllib.request.build_opener()
        opener.addheaders = [("User-agent", "Mozilla/5.0")]
        urllib.request.install_opener(opener)
        response = urllib.request.urlretrieve(api_url)
    except KeyboardInterrupt:
        # when CTRL+C is pressed during the execution of this script
        logging.error(user_interrupt_text)
        raise ValueError(user_interrupt_text)

    # Make the base directories for this GitHub folder
    (output_dir / download_dirs).mkdir(parents=True, exist_ok=True)

    with open(response[0], "r") as f:
        data = json.load(f)

    for entry in data:
        # if the entry is a further folder, walk through it
        if entry["type"] == "dir":
            download_github_folder(repo_url=entry["html_url"],
                                   output_dir=output_dir)

        # if the entry is a file, download it
        elif entry["type"] == "file":
            try:
                opener = urllib.request.build_opener()
                opener.addheaders = [("User-agent", "Mozilla/5.0")]
                urllib.request.install_opener(opener)
                # download the file
                save_path = output_dir / entry["path"]
                urllib.request.urlretrieve(entry["download_url"], str(save_path))
                logging.info("Downloaded: %s", entry["path"])

            except KeyboardInterrupt:
                # when CTRL+C is pressed during the execution of this script
                logging.error(user_interrupt_text)
                raise ValueError(user_interrupt_text)
