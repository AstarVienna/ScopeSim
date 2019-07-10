import shutil
import os
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
        rc keyword FILE_SERVER_BASE_URL

    Returns
    -------
    all_pkgs : list of str
        A list of paths to the ``.zip`` packages relative to ``url``

    """
    if url is None:
        url = rc.__rc__["FILE_SERVER_BASE_URL"]

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
        If left as None, defaults to the rc keyword FILE_LOCAL_DOWNLOADS_PATH

    url : str
        The URL of the IRDB HTTP server. If left as None, defaults to the
        rc keyword FILE_SERVER_BASE_URL

    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

    """
    if url is None:
        url = rc.__rc__["FILE_SERVER_BASE_URL"]

    if save_dir is None:
        save_dir = rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"]

    try:
        file_path = download_file(url + pkg_path, cache=True)
        save_path = save_dir + os.path.basename(pkg_path)
        shutil.copy2(file_path, save_path)
    except HTTPError:
        ValueError("Unable to find file: {}".format(url + pkg_path))

    return os.path.abspath(save_path)
