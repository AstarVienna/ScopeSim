# -*- coding: utf-8 -*-
"""
Used only by the `database` submodule.

Original comment for these functions:
    2022-04-10 (KL)
    Code taken directly from https://github.com/sdushantha/gitdir
    Adapted for ScopeSim usage.
    Many thanks to the authors!

"""

import logging
import re
from pathlib import Path
from typing import Union

import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

from .download_utils import initiate_download, handle_download


HTTP_RETRY_CODES = [403, 404, 429, 500, 501, 502, 503]


class ServerError(Exception):
    """Some error with the server or connection to the server."""


def create_github_url(url: str) -> None:
    """
    From the given url, produce a URL that is compatible with Github's REST API.

    Can handle blob or tree paths.
    """
    repo_only_url = re.compile(r"https:\/\/github\.com\/[a-z\d](?:[a-z\d]|-(?=[a-z\d])){0,38}\/[a-zA-Z0-9]+$")
    re_branch = re.compile("/(tree|blob)/(.+?)/")

    # Check if the given url is a url to a GitHub repo. If it is, tell the
    # user to use 'git clone' to download it
    if re.match(repo_only_url,url):
        message = ("✘ The given url is a complete repository. Use 'git clone'"
                   " to download the repository")
        logging.error(message)
        raise ValueError(message)

    # extract the branch name from the given url (e.g master)
    branch = re_branch.search(url)
    download_dirs = url[branch.end():]
    api_url = (url[:branch.start()].replace("github.com", "api.github.com/repos", 1) +
               f"/contents/{download_dirs}?ref={branch.group(2)}")
    return api_url, download_dirs


def download_github_folder(repo_url: str,
                           output_dir: Union[Path, str] = "./") -> None:
    """
    Downloads the files and directories in repo_url.

    Re-written based on the on the download function
    `here <https://github.com/sdushantha/gitdir/blob/f47ce9d85ee29f8612ce5ae804560a12b803ddf3/gitdir/gitdir.py#L55>`_
    """
    output_dir = Path(output_dir)

    # convert repo_url into an api_url
    api_url, download_dirs = create_github_url(repo_url)

    # get the contents of the github folder
    try:
        retry_strategy = Retry(total=3, backoff_factor=2,
                               status_forcelist=HTTP_RETRY_CODES,
                               allowed_methods=["GET"])
        adapter = HTTPAdapter(max_retries=retry_strategy)
        with requests.Session() as session:
            session.mount("https://", adapter)
            data = session.get(api_url).json()
    except (requests.exceptions.ConnectionError,
            requests.exceptions.RetryError) as error:
        logging.error(error)
        raise ServerError("Cannot connect to server. "
                          f"Attempted URL was: {api_url}.") from error
    except Exception as error:
        logging.error(("Unhandled exception occured while accessing server."
                      "Attempted URL was: %s."), api_url)
        logging.error(error)
        raise error

    # Make the base directories for this GitHub folder
    (output_dir / download_dirs).mkdir(parents=True, exist_ok=True)

    for entry in data:
        # if the entry is a further folder, walk through it
        if entry["type"] == "dir":
            download_github_folder(repo_url=entry["html_url"],
                                   output_dir=output_dir)

        # if the entry is a file, download it
        elif entry["type"] == "file":
            try:
                # download the file
                save_path = output_dir / entry["path"]
                response = initiate_download(entry["download_url"])
                handle_download(response, save_path, entry["path"],
                                padlen=0, disable_bar=True)
                logging.info("Downloaded: %s", entry["path"])

            except (requests.exceptions.ConnectionError,
                    requests.exceptions.RetryError) as error:
                logging.error(error)
                raise ServerError("Cannot connect to server. "
                                  f"Attempted URL was: {api_url}.") from error
            except Exception as error:
                logging.error(("Unhandled exception occured while accessing "
                              "server. Attempted URL was: %s."), api_url)
                logging.error(error)
                raise error
