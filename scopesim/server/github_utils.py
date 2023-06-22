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
import json
from pathlib import Path
from typing import Union

import urllib

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
    user_interrupt_text = "GitHub download interrupted by User"
    try:
        opener = urllib.request.build_opener()
        opener.addheaders = [("User-agent", "Mozilla/5.0")]
        urllib.request.install_opener(opener)
        response = urllib.request.urlretrieve(api_url)
    except KeyboardInterrupt as error:
        logging.error(user_interrupt_text)
        raise error

    # Make the base directories for this GitHub folder
    (output_dir / download_dirs).mkdir(parents=True, exist_ok=True)

    with open(response[0], "r") as f:
        data = json.load(f)

    for entry in data:
        # if the entry is a further folder, walk through it
        if entry["type"] == "dir":
            download_github_folder(repo_url=entry["html_url"],
                                   output_dir=output_dir)

        # if the entry is a file, download it
        elif entry["type"] == "file":
            try:
                opener = urllib.request.build_opener()
                opener.addheaders = [("User-agent", "Mozilla/5.0")]
                urllib.request.install_opener(opener)
                # download the file
                save_path = output_dir / entry["path"]
                urllib.request.urlretrieve(entry["download_url"], str(save_path))
                logging.info("Downloaded: %s", entry["path"])

            except KeyboardInterrupt as error:
                logging.error(user_interrupt_text)
                raise error
