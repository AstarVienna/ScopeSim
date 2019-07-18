import shutil
import os
import zipfile
from urllib3.exceptions import HTTPError

import requests
import bs4
from astropy.utils.data import download_file

from scopesim import rc


def get_server_elements(url, unique_str="/"):
    """
    Returns a list of file and/or directory paths on the HTTP server ``url``

    Parameters
    ----------
    url : str
        The URL of the IRDB HTTP server.

    unique_str : str
        A unique string to look for in the beautiful HTML soup:
        "/" for directories this, ".zip" for packages

    Returns
    -------
    paths : list
        List of paths containing in ``url`` which contain ``unique_str``

    """
    try:
        result = requests.get(url).content
    except:
        raise ValueError("URL returned error: {}".format(url))

    soup = bs4.BeautifulSoup(result)
    paths = soup.findAll("a", href=True)
    paths = [tmp.string for tmp in paths if unique_str in tmp.string]

    return paths


def list_packages(url=None):
    """
    List all ``.zip`` packages found under ``url``

    Parameters
    ----------
    url : str
        The URL of the IRDB HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    Returns
    -------
    all_pkgs : list of str
        A list of paths to the ``.zip`` packages relative to ``url``

    """
    if url is None:
        url = rc.__config__["!SIM.file.server_base_url"]

    all_pkgs = []
    folders = get_server_elements(url, "/")
    for folder in folders:
        pkgs = get_server_elements(url + folder, ".zip")
        all_pkgs += [folder + pkg for pkg in pkgs]

    return all_pkgs


def download_package(pkg_path, save_dir=None, url=None):
    """
    Downloads a package to the local disk

    Parameters
    ----------
    pkg_path : str
        A ``.zip`` package path as given by ``list_packages()``

    save_dir : str
        The place on the local disk where the ``.zip`` package is to be saved.
        If left as None, defaults to the value in
        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    url : str
        The URL of the IRDB HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

    """
    if url is None:
        url = rc.__config__["!SIM.file.server_base_url"]
    if save_dir is None:
        save_dir = rc.__config__["!SIM.file.local_packages_path"]

    try:
        use_cached_file = rc.__config__["!SIM.file.use_cached_downloads"]
        cache_path = download_file(url + pkg_path, cache=use_cached_file)
        save_path = os.path.join(save_dir, os.path.basename(pkg_path))
        file_path = shutil.copy2(cache_path, save_path)

        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(save_dir)

    except HTTPError:
        ValueError("Unable to find file: {}".format(url + pkg_path))

    return os.path.abspath(save_path)
