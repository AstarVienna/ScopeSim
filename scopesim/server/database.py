"""
Functions to download instrument packages and example data
"""
import shutil
import os
import zipfile
import logging
from urllib3.exceptions import HTTPError

import yaml
import requests
import bs4
from astropy.utils.data import download_file

from scopesim import rc
from .gitdir import download as download_github_folder


def get_server_package_list():
    url = rc.__config__["!SIM.file.server_base_url"]
    response = requests.get(url + "packages.yaml")
    pkgs_dict = yaml.full_load(response.text)

    return pkgs_dict


def get_server_folder_contents(dir_name, unique_str=".zip"):
    url = rc.__config__["!SIM.file.server_base_url"] + dir_name

    try:
        result = requests.get(url).content
    except Exception as error:
        raise ValueError(f"URL returned error: {url}") from error

    soup = bs4.BeautifulSoup(result, features="lxml")
    hrefs = soup.findAll("a", href=True)
    pkgs = [href.string for href in hrefs
            if href.string is not None and ".zip" in href.string]

    return pkgs


def list_packages(pkg_name=None):
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
    pkgs_dict = get_server_package_list()

    if pkg_name is None:
        pkg_names = list(pkgs_dict.keys())
    elif pkg_name in pkgs_dict:
        path = pkgs_dict[pkg_name]["path"]
        pkgs = get_server_folder_contents(path)
        pkg_names = [pkg for pkg in pkgs if pkg_name in pkg]

    return pkg_names


def download_packages(pkg_names, release="stable", save_dir=None, from_cache=None):
    """
    Download one or more packages to the local disk

    1. Download stable, dev
    2. Download specific version
    3. Download from github via url

    Parameters
    ----------
    pkg_names : str, list
        A list of package name, see ``list_packages()``

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

        # Specific package from a Gtihub commit hash (use "@" or ":")
        download_packages("ELT", release="github:728761fc76adb548696205139e4e9a4260401dfc")
        download_packages("ELT", release="github@728761fc76adb548696205139e4e9a4260401dfc")

    """
    base_url = rc.__config__["!SIM.file.server_base_url"]

    pkgs_dict = get_server_package_list()

    if isinstance(pkg_names, str):
        pkg_names = [pkg_names]

    save_paths = []
    for pkg_name in pkg_names:
        if pkg_name in pkgs_dict:
            pkg_dict = pkgs_dict[pkg_name]
            path = pkg_dict["path"] + "/"

            from_github = False
            if release in ["stable", "latest"]:
                zip_name = pkg_dict[release]
                pkg_url = f"{base_url}{path}/{zip_name}.zip"
            elif "github" in release:
                base_url = "https://github.com/AstarVienna/irdb/tree/"
                github_hash = release.split(":")[-1].split("@")[-1]
                pkg_url = f"{base_url}{github_hash}/{pkg_name}"
                from_github = True
            else:
                zip_name = f"{pkg_name}.{release}.zip"
                pkg_variants = get_server_folder_contents(path)
                if zip_name not in pkg_variants:
                    raise ValueError(f"{zip_name} is not amoung the hosted "
                                     f"variants: {pkg_variants}")
                pkg_url = f"{base_url}{path}/{zip_name}"

            if save_dir is None:
                save_dir = rc.__config__["!SIM.file.local_packages_path"]
            if not os.path.exists(save_dir):
                os.mkdir(save_dir)

            if not from_github:
                try:
                    if from_cache is None:
                        from_cache = rc.__config__["!SIM.file.use_cached_downloads"]
                    cache_path = download_file(pkg_url, cache=from_cache)
                    save_path = os.path.join(save_dir, f"{pkg_name}.zip")
                    file_path = shutil.copy2(cache_path, save_path)

                    with zipfile.ZipFile(file_path, 'r') as zip_ref:
                        zip_ref.extractall(save_dir)

                except HTTPError:
                    ValueError(f"Unable to find file: {url + pkg_path}")
            else:
                download_github_folder(repo_url=pkg_url, output_dir=save_dir)
                save_path = save_dir

            save_paths += [os.path.abspath(save_path)]

    return save_paths

# ==============================================================================
# Funtions below from from OLD_database.py
# ==============================================================================

# for backwards compatibility
def download_package(pkg_path, save_dir=None, url=None, from_cache=None):
    """
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
    # todo: add proper depreciation warning
    text = "Function Depreciated --> please use scopesim.download_package-s-()"
    logging.warning(text)
    print(text)

    pkg_names = [pkg.split("/")[-1].replace(".zip", "") for pkg in pkg_path]
    return download_packages(pkg_names, release="stable", save_dir=save_dir,
                             from_cache=from_cache)

def get_server_elements(url, unique_str="/"):
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


def list_example_data(url=None, return_files=False, silent=False):
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


def download_example_data(file_path, save_dir=None, url=None, from_cache=None):
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
    save_path : str
        The absolute path to the saved files
    """
    if isinstance(file_path, (list, tuple)):
        save_path = [download_example_data(thefile, save_dir, url)
                     for thefile in file_path]
    elif isinstance(file_path, str):

        if url is None:
            url = rc.__config__["!SIM.file.server_base_url"]
        if save_dir is None:
            save_dir = os.getcwd()
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        try:
            if from_cache is None:
                from_cache = rc.__config__["!SIM.file.use_cached_downloads"]
            cache_path = download_file(url + "example_data/" + file_path,
                                       cache=from_cache)
            save_path = os.path.join(save_dir, os.path.basename(file_path))
            file_path = shutil.copy2(cache_path, save_path)
        except HTTPError:
            ValueError(f"Unable to find file: {url + 'example_data/' + file_path}")

        save_path = os.path.abspath(save_path)

    return save_path