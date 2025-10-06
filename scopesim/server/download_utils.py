# -*- coding: utf-8 -*-
"""Used only by the `database` and `github_utils` submodules."""

import re
from collections.abc import Iterator
from zipfile import ZipFile
from pathlib import Path

import httpx
import bs4
import pooch

from tqdm.auto import tqdm
# from tqdm.contrib.logging import logging_redirect_tqdm
# put with logging_redirect_tqdm(loggers=all_loggers): around tqdm
# Note: seems to work without that so far...

from .. rc import __config__
from ..utils import get_logger


logger = get_logger(__name__)


class ServerError(Exception):
    """Some error with the server or connection to the server."""


def _get_svr_conf():
    return __config__["!SIM.file"]


def get_base_url():
    """Get instrument package server URL from rc.__config__."""
    return _get_svr_conf()["server_base_url"]


def get_local_packages_path():
    """Get local instrument package save directory from rc.__config__."""
    return Path(_get_svr_conf()["local_packages_path"])


def _make_tqdm_kwargs(desc: str = ""):
    # width, _ = get_terminal_size((50, 20))
    # bar_width = max(int(.8 * width) - 30 - len(desc), 10)
    tqdm_kwargs = {
        # "bar_format": f"{{l_bar}}{{bar:{bar_width}}}{{r_bar}}{{bar:-{bar_width}b}}",
        "colour": "green",
        "desc": desc
        }
    return tqdm_kwargs


def create_client(
    base_url: str,
    cached: bool = False,
    cache_name: str = "",
) -> httpx.Client:
    """
    Create httpx Client instance, should support cache at some point.

    Parameters
    ----------
    base_url : str
        Server base URL.
    cached : bool, optional
        Not yet implemented. The default is False.
    cache_name : str, optional
        Not yet implemented. The default is "".

    Returns
    -------
    client : httpx.Client
        Client instance.

    """
    if cached:
        raise NotImplementedError("Caching not yet implemented with httpx.")
    transport = httpx.HTTPTransport(retries=5)
    client = httpx.Client(base_url=base_url, timeout=2, transport=transport)
    return client


def handle_download(
    client: httpx.Client,
    sub_url: str,
    save_path: Path,
    name: str = "",
    padlen: int = 0,
    chunk_size: int = 128,
    disable_bar: bool = False,
    params: dict | None = None,
) -> None:
    """
    Perform a streamed download and write the content to disk.

    Parameters
    ----------
    client : httpx.Client
        Client instance.
    sub_url : str
        URL of file to be downloaded, relative to the `client`'s base URL.
    save_path : Path
        Path including file name to save the downloaded file to.
    name : str, optional
        Optional display name for progress bar. Has no influence on saved
        filename. The default is "".
    padlen : int, optional
        Optional padding for `name` in progress bar. Useful to get a nicer
        output when doenloading multiple files. The default is 0.
    chunk_size : int, optional
        Chunk size for download stream. The default is 128.
    disable_bar : bool, optional
        If True, no progress bar is printed. The default is False.
    params : dict | None, optional
        Any query parameters to include in the URL, will be forwarded to the
        ``send_get()`` call. The default is None.

    Raises
    ------
    ServerError
        Raised if any error occurs during the process.

    Returns
    -------
    None

    """
    tqdm_kwargs = _make_tqdm_kwargs(f"Downloading {name:<{padlen}}")
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    stream = send_get(client, sub_url, True, params)

    try:
        with stream as response:
            response.raise_for_status()
            total = int(response.headers.get("Content-Length", 0))

            with (save_path.open("wb") as file_outer,
                  tqdm.wrapattr(
                    file_outer, "write", miniters=1, total=total,
                    **tqdm_kwargs, disable=disable_bar) as file_inner):
                for chunk in response.iter_bytes(chunk_size=chunk_size):
                    file_inner.write(chunk)

    except httpx.HTTPStatusError as err:
        logger.error("Error response %s while requesting %s.",
                     err.response.status_code, err.request.url)
        raise ServerError("Cannot connect to server.") from err
    except OSError as err:
        logger.error("Error %s attempting to access path %s.", err, save_path)
        raise ServerError("Cannot access save path.") from err
    except ServerError:
        raise
    except Exception as err:
        logger.exception("Unhandled exception while accessing server.")
        raise ServerError("Cannot perform download.") from err


def handle_unzipping(save_path: Path, save_dir: Path,
                     pkg_name: str, padlen: int) -> None:
    """Unpack a zipped folder, usually called right after downloading."""
    with ZipFile(save_path, "r") as zip_ref:
        namelist = zip_ref.namelist()
        tqdm_kwargs = _make_tqdm_kwargs(f"Extracting  {pkg_name:<{padlen}}")
        for file in tqdm(iterable=namelist, total=len(namelist), **tqdm_kwargs):
            zip_ref.extract(file, save_dir)


def send_get(
    client: httpx.Client,
    sub_url: str,
    stream: bool = False,
    params: dict | None = None,
) -> httpx.Response:
    """
    Send a GET request (streamed or not) using an existing client.

    The point of this function is mostly elaborate exception handling.

    Parameters
    ----------
    client : httpx.Client
        Client instance.
    sub_url : str
        Sub-URL to be appended to `client.base_url`.
    stream : bool, optional
        Whether to stream the response. The default is False.
    params : dict | None, optional
        Any query parameters to include in the URL. The default is None.

    Raises
    ------
    ServerError
        Raised if any error occurs during the process.

    Returns
    -------
    response : httpx.Response
        Server response.

    """
    try:
        if stream:
            response = client.stream("GET", sub_url, params=params)
        else:
            response = client.get(sub_url, params=params)
            response.raise_for_status()
    except httpx.RequestError as err:
        logger.exception("An error occurred while requesting %s.",
                         err.request.url)
        raise ServerError("Cannot connect to server.") from err
    except httpx.HTTPStatusError as err:
        logger.error("Error response %s while requesting %s.",
                     err.response.status_code, err.request.url)
        raise ServerError("Cannot connect to server.") from err
    except Exception as err:
        logger.exception("Unhandled exception while accessing server.")
        raise ServerError("Cannot connect to server.") from err

    return response


def get_server_folder_contents(client, dir_name: str,
                               unique_str: str = ".zip$") -> Iterator[str]:
    """Find all zip files in a given server folder."""
    dir_name = dir_name + "/" if not dir_name.endswith("/") else dir_name
    response = send_get(client, dir_name)

    soup = bs4.BeautifulSoup(response.content, features="lxml")
    hrefs = soup.find_all("a", href=True, string=re.compile(unique_str))
    pkgs = (href.string for href in hrefs)

    return pkgs


def create_retriever(
    collection: str,
    save_dir: Path | str | None = None,
    url: str | None = None,
) -> pooch.Pooch:
    """
    Create Pooch retriever and load example data registry.

    This functions is meant for internal use. Currently supported values for
    `collection` are "example_data" (used by the example data download
    functions, which *are* part of the user-facing API), "atmo" (for
    atmospheric libraries) and "psfs" (not yet implemented). The main purpose
    of this download infrastructure is to be able to store (relatively) large
    static files in a remote location, instead of having to bundle them in the
    package (or in IRDB packages).

    The ``Pooch`` instance returned by this function can then download any file
    listed in the respective registry via ``retriever.fetch(<filename>,
    progressbar=True)``, which returns the full path to the locally cached file.
    Use ``retriever.registry_files`` to display which files are available for
    download in each instance.

    Parameters
    ----------
    collection : str
        "Server folder" containing data to be retrieved.
    save_dir : Path | str | None, optional
        Local cache location. If None (the default), use the one specified in
        the global configuration file.
    url : str | None, optional
        Server base URL. If None (the default), use the one specified in
        the global configuration file.

    Returns
    -------
    retriever : pooch.Pooch
        Pooch retriever instance.

    """
    svrconf = _get_svr_conf()

    url = url or (get_base_url() + svrconf[collection]["suburl"])
    save_dir = save_dir or (Path.home() / ".astar/scopesim" / collection)
    save_dir = Path(save_dir)

    retriever = pooch.create(
        path=save_dir,
        base_url=url,
        retry_if_failed=svrconf["retries"],
    )
    registry_file = Path(__file__).parent / svrconf[collection]["hash_file"]
    retriever.load_registry(registry_file)

    return retriever
