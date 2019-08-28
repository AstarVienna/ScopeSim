import shutil
import os
import zipfile
from urllib3.exceptions import HTTPError
import glob

import requests
import bs4
from astropy.utils.data import download_file

from scopesim import rc


def get_local_packages(path):
    """
    List the packages that are available in the directory ``path``

    Parameters
    ----------
    path : str
        Directory containing all local instrument package files

    Returns
    -------
    pkgs : list
        Names of packages on the local disk

    """
    dirnames = os.listdir(path)
    pkgs = []

    for dname in dirnames:
        if os.path.exists(os.path.join(path, dname, dname+".yaml")):
            pkgs += [dname]

    return pkgs


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


def list_packages(location="all", url=None, local_dir=None,
                  return_pkgs=False, silent=False):
    """
    List all ``.zip`` packages found under ``url``

    Parameters
    ----------
    location : str
        ["server", "local", "all"] To look for packages on the server or on the
        local hard-drive, or both

    url : str
        The URL of the IRDB HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    local_dir : str
        Path to the folder where the local packages are stored (or downloaded)
        If left as None, defaults to
        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    return_pkgs : bool
        If True, returns a list of package names

    silent : bool
        If True, does not print the list of package names

    Returns
    -------
    all_pkgs : list of str
        A list of paths to the ``.zip`` packages relative to ``url``
        The full string should be passed to download_package

    """
    if url is None:
        url = rc.__config__["!SIM.file.server_base_url"]
    if local_dir is None:
        local_dir = rc.__config__["!SIM.file.local_packages_path"]

    server_pkgs = []
    folders = get_server_elements(url, "/")
    for folder in folders:
        pkgs = get_server_elements(url + folder, ".zip")
        server_pkgs += [folder + pkg for pkg in pkgs]

    local_pkgs = get_local_packages(local_dir)

    def print_package_list(the_pkgs, loc=""):
        print("\nPackages saved {}".format(loc) + "\n" + "=" * (len(loc) + 15))
        for pkg in the_pkgs:
            print(pkg)

    if "local" in location:
        if not silent:
            print_package_list(local_pkgs, "locally: {}".format(local_dir))
        if return_pkgs:
            return local_pkgs

    elif "server" in location:
        if not silent:
            print_package_list(server_pkgs, "on the server: {}".format(url))
        if return_pkgs:
            return local_pkgs

    elif "all" in location:
        if not silent:
            print_package_list(server_pkgs, "on the server: {}".format(url))
            print_package_list(local_pkgs, "locally: {}".format(local_dir))
        if return_pkgs:
            return server_pkgs, local_pkgs


def download_package(pkg_path, save_dir=None, url=None):
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

    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

    """
    if isinstance(pkg_path, (list, tuple)):
        save_path = [download_package(pkg, save_dir, url) for pkg in pkg_path]

    elif isinstance(pkg_path, str):
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

        save_path = os.path.abspath(save_path)

    return save_path
