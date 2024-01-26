# -*- coding: utf-8 -*-
"""
Used only by the `database` submodule.

Original comment for these functions:
    2022-04-10 (KL)
    Code taken directly from https://github.com/sdushantha/gitdir
    Adapted for ScopeSim usage.
    Many thanks to the authors!

"""

import re
from pathlib import Path
from typing import Union

from .download_utils import handle_download, send_get, create_client
from ..utils import get_logger


logger = get_logger(__name__)


def create_github_url(url: str) -> None:
    """
    From the given url, produce a URL compatible with Github's REST API.

    Can handle blob or tree paths.
    """
    repo_only_url = re.compile(r"https:\/\/github\.com\/[a-z\d](?:[a-z\d]|-(?=[a-z\d])){0,38}\/[a-zA-Z0-9]+$")
    re_branch = re.compile("/(tree|blob)/(.+?)/")

    # Check if the given url is a url to a GitHub repo. If it is, tell the
    # user to use 'git clone' to download it
    if re.match(repo_only_url, url):
        message = ("✘ The given url is a complete repository. Use 'git clone'"
                   " to download the repository")
        logger.error(message)
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
    Download the files and directories in repo_url.

    Re-written based on the on the download function
    `here <https://github.com/sdushantha/gitdir/blob/f47ce9d85ee29f8612ce5ae804560a12b803ddf3/gitdir/gitdir.py#L55>`_
    """
    output_dir = Path(output_dir)

    # convert repo_url into an api_url
    api_url, download_dirs = create_github_url(repo_url)

    # get the contents of the github folder
    with create_client("", cached=False) as client:
        data = send_get(client, api_url).json()

        # Make the base directories for this GitHub folder
        (output_dir / download_dirs).mkdir(parents=True, exist_ok=True)

        for entry in data:
            # if the entry is a further folder, walk through it
            if entry["type"] == "dir":
                download_github_folder(repo_url=entry["html_url"],
                                       output_dir=output_dir)

            # if the entry is a file, download it
            elif entry["type"] == "file":
                # download the file
                save_path = output_dir / entry["path"]
                handle_download(client, entry["download_url"], save_path,
                                entry["path"], padlen=0, disable_bar=True)
                logger.info("Downloaded: %s", entry["path"])
