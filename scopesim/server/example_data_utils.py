# -*- coding: utf-8 -*-
"""
Store the example data functions here instead of polluting database.py
"""

import shutil
from pathlib import Path
from typing import List, Optional, Union, Iterable

from urllib.error import HTTPError
from urllib3.exceptions import HTTPError as HTTPError3

import requests
import bs4

from astropy.utils.data import download_file

from scopesim import rc

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
