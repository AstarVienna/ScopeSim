# -*- coding: utf-8 -*-
"""
Used only by the `database` and `github_utils` submodules.
"""

from zipfile import ZipFile
from pathlib import Path
from shutil import get_terminal_size

import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from requests_cache import CachedSession
from tqdm import tqdm
# from tqdm.contrib.logging import logging_redirect_tqdm
# put with logging_redirect_tqdm(loggers=all_loggers): around tqdm



HTTP_RETRY_CODES = [403, 404, 429, 500, 501, 502, 503]


def _make_tqdm_kwargs(desc: str = ""):
    width, _ = get_terminal_size((50, 20))
    bar_width = max(int(.8 * width) - 30 - len(desc), 10)
    tqdm_kwargs = {
        "bar_format": f"{{l_bar}}{{bar:{bar_width}}}{{r_bar}}{{bar:-{bar_width}b}}",
        "colour": "green",
        "desc": desc
        }
    return tqdm_kwargs


def _create_session(cached: bool = False, cache_name: str = ""):
    if cached:
        return CachedSession(cache_name)
    return requests.Session()


def initiate_download(pkg_url: str,
                      cached: bool = False, cache_name: str = "",
                      total: int = 5, backoff_factor: int = 2):
    retry_strategy = Retry(total=total, backoff_factor=backoff_factor,
                           status_forcelist=HTTP_RETRY_CODES,
                           allowed_methods=["GET"])
    adapter = HTTPAdapter(max_retries=retry_strategy)
    with _create_session(cached, cache_name) as session:
        session.mount("https://", adapter)
        response = session.get(pkg_url, stream=True)
    return response


def handle_download(response, save_path: Path, pkg_name: str,
                    padlen: int, chunk_size: int = 128,
                    disable_bar=False) -> None:
    tqdm_kwargs = _make_tqdm_kwargs(f"Downloading {pkg_name:<{padlen}}")
    total = int(response.headers.get("content-length", 0))
    # Turn this into non-nested double with block in Python 3.9 or 10 (?)
    with save_path.open("wb") as file_outer:
        with tqdm.wrapattr(file_outer, "write", miniters=1, total=total,
                           **tqdm_kwargs, disable=disable_bar) as file_inner:
            for chunk in response.iter_content(chunk_size=chunk_size):
                file_inner.write(chunk)


def handle_unzipping(save_path: Path, save_dir: Path,
                     pkg_name: str, padlen: int) -> None:
    with ZipFile(save_path, "r") as zip_ref:
        namelist = zip_ref.namelist()
        tqdm_kwargs = _make_tqdm_kwargs(f"Extracting  {pkg_name:<{padlen}}")
        for file in tqdm(iterable=namelist, total=len(namelist), **tqdm_kwargs):
            zip_ref.extract(file, save_dir)
