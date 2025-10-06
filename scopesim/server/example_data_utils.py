# -*- coding: utf-8 -*-
"""Convenienve functions for listing and downloading example datasets."""

from pathlib import Path

from .download_utils import create_retriever


def list_example_data(
    url: str | None = None,
    return_files: bool = False,
    silent: bool = False,
) -> list[str]:
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
    retriever = create_retriever("example_data", url=url)
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


def download_example_data(
    *files: str,
    save_dir: Path | str | None = None,
    url: str | None = None,
) -> list[Path]:
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

    .. versionchanged:: 0.8.4

       Passing a list to ``download_example_data`` is deprecated since version
       0.8.4, this function now accepts multiple file names in *args-style.

    .. versionchanged:: 0.11.0

       Passing a list to ``download_example_data`` as the first argument will
       now throw a TypeError. This is to catch any remaining uses of the old
       call signature of this function. From version 0.12 onwards, agruments
       will be silently passed to the `retriever`.

    """
    if isinstance(files[0], list):
        raise TypeError(
            "Passing a list to download_example_data is deprecated. "
            "Simply pass filenames as *args, i.e. "
            "download_example_data(\"foo.fits\", \"bar.fits\")."
        )

    retriever = create_retriever("example_data", save_dir, url)

    save_paths = []
    for fname in files:
        save_paths.append(Path(retriever.fetch(fname, progressbar=True)))

    return save_paths
