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


    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

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

            if release in ["stable", "latest"]:
                zip_name = pkg_dict[release]
                pkg_url = f"{base_url}{path}/{zip_name}.zip"
            elif "github" in release:
                # try to get from GitHub
                pass
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

            save_paths += [os.path.abspath(save_path)]

    return save_paths