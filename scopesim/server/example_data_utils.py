# -*- coding: utf-8 -*-
"""Store the example data functions here instead of polluting database.py."""

from warnings import warn
from pathlib import Path
from typing import Optional, Union

import pooch

from scopesim import rc


def _create_retriever(save_dir: Optional[Union[Path, str]] = None,
                      url: Optional[str] = None) -> pooch.Pooch:
    """Create Pooch retriever and load example data registry."""
    svrconf = rc.__config__["!SIM.file"]

    url = url or (svrconf["server_base_url"] + svrconf["example_data_suburl"])
    save_dir = save_dir or (Path.home() / ".astar/scopesim")
    save_dir = Path(save_dir)

    retriever = pooch.create(
        path=save_dir,
        base_url=url,
        retry_if_failed=3)
    registry_file = Path(__file__).parent / svrconf["example_data_hash_file"]
    retriever.load_registry(registry_file)

    return retriever


def list_example_data(url: Optional[str] = None,
                      return_files: bool = False,
                      silent: bool = False) -> list[str]:
    """
    List all example files found under ``url``.

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
    retriever = _create_retriever(url=url)
    server_files = retriever.registry_files

    if not silent:
        linelen = max(len(max(server_files, key=len)),
                      len(retriever.base_url) + 6, 36)
        print(f"\n{f' Files available on the server: ':=^{linelen}}")
        print(f"{f' {retriever.base_url} ':=^{linelen}}\n")
        for file in server_files:
            print(file)

    if return_files:
        return server_files

    return None


def download_example_data(*files: str,
                          save_dir: Optional[Union[Path, str]] = None,
                          url: Optional[str] = None) -> list[Path]:
    """
    Download example fits files to the local disk.

    Parameters
    ----------
    files : str(s)
        Name(s) of FITS file(s) as given by ``list_example_data()``.

    save_dir : str
        The place on the local disk where the downloaded files are to be saved.
        If left as None, defaults to the "~/.astar/scopesim".

    url : str
        The URL of the database HTTP server. If left as None, defaults to the
        value in scopesim.rc.__config__["!SIM.file.server_base_url"]

    Returns
    -------
    save_path : list of Paths
        The absolute path(s) to the saved files
    """
    if isinstance(files[0], list):  # to preserve combatibility with notebooks
        warn("Passing a list to download_example_data is deprecated. "
             "Simply pass filenames as *args, i.e. "
             "download_example_data(\"foo.fits\", \"bar.fits\").",
             DeprecationWarning, stacklevel=2)
        files = files[0]

    retriever = _create_retriever(save_dir, url)

    save_paths = []
    for fname in files:
        save_paths.append(Path(retriever.fetch(fname, progressbar=True)))

    return save_paths
