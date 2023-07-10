"""
Functions to download instrument packages and example data
"""
import shutil
import os
import zipfile
import logging

from scopesim.server.database import get_server_elements
from urllib3.exceptions import HTTPError

from astropy.utils.data import download_file

from scopesim import rc

from warnings import warn

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
    warn("Function Depreciated --> please use scopesim.download_package-s-()",
         DeprecationWarning, stacklevel=2)
    dirnames = os.listdir(path)
    pkgs = []

    for dname in dirnames:
        if os.path.exists(os.path.join(path, dname, dname+".yaml")):
            pkgs += [dname]

    return pkgs


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

    def print_package_list(the_pkgs, loc=""):
        print(f"\nPackages saved {loc}\n" + "=" * (len(loc) + 15))
        for pkg in the_pkgs:
            print(pkg)

    if url is None:
        url = rc.__config__["!SIM.file.server_base_url"]
    if local_dir is None:
        local_dir = rc.__config__["!SIM.file.local_packages_path"]

    return_pkgs_list = []

    if location.lower() in ["local", "all"]:
        local_pkgs = get_local_packages(local_dir)
        if not silent:
            print_package_list(local_pkgs, f"locally: {local_dir}")
            return_pkgs_list += local_pkgs

    if location.lower() in ["server", "all"]:
        server_pkgs = []
        folders = get_server_elements(url, "/")
        for folder in folders:
            pkgs = get_server_elements(url + folder, ".zip")
            server_pkgs += [folder + pkg for pkg in pkgs]
        if not silent:
            print_package_list(server_pkgs, f"on the server: {url}")
            return_pkgs_list += server_pkgs

    if return_pkgs:
        return return_pkgs_list

    return None

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
            raise ValueError(f"Unable to find file: {url + pkg_path}")

        save_path = os.path.abspath(save_path)

    return save_path
