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


def download_packages(pkg_names, save_dir=None, url=None, from_cache=None):
    """
    Download one or more packages to the local disk

    Parameters
    ----------
    pkg_names : str, list
        A list of package name, see ``list_packages()``

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
    if isinstance(pkg_path, (list, tuple)):
        save_path = [download_package(pkg, save_dir, url) for pkg in pkg_path]

    elif isinstance(pkg_path, str):
        if pkg_path[-4:] != ".zip":
            logging.warning("Appended '.zip' to %s", pkg_path)
            pkg_path += ".zip"

        if url is None:
            url = rc.__config__["!SIM.file.server_base_url"]
        if save_dir is None:
            save_dir = rc.__config__["!SIM.file.local_packages_path"]

        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        try:
            if from_cache is None:
                from_cache = rc.__config__["!SIM.file.use_cached_downloads"]
            cache_path = download_file(url + pkg_path, cache=from_cache)
            save_path = os.path.join(save_dir, os.path.basename(pkg_path))
            file_path = shutil.copy2(cache_path, save_path)

            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(save_dir)

        except HTTPError:
            ValueError(f"Unable to find file: {url + pkg_path}")

        save_path = os.path.abspath(save_path)

    return save_path


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
    pkgs_dict = get_server_package_list()

    if pkg_name is None:
        pkg_names = list(pkgs_dict.keys())
    elif pkg_name in pkgs_dict:
        path = pkgs_dict[pkg_name]["path"]
        pkgs = get_server_folder_contents(path)
        pkg_names = [pkg for pkg in pkgs if pkg_name in pkg]

    return pkg_names
