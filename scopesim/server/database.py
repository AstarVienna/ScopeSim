"""
Functions to download instrument packages and example data
"""
import json
import re
import shutil
import os
import urllib.request
import zipfile
import logging
from urllib3.exceptions import HTTPError

import yaml
import requests
import bs4
from astropy.utils.data import download_file

from scopesim import rc


def get_server_package_list():
    url = rc.__config__["!SIM.file.server_base_url"]
    response = requests.get(url + "packages.yaml")
    pkgs_dict = yaml.full_load(response.text)

    return pkgs_dict


def get_server_folder_contents(dir_name, unique_str=".zip"):
    url = rc.__config__["!SIM.file.server_base_url"] + dir_name

    try:
        result = requests.get(url).content
    except Exception as error:
        raise ValueError(f"URL returned error: {url}") from error

    soup = bs4.BeautifulSoup(result, features="lxml")
    hrefs = soup.findAll("a", href=True)
    pkgs = [href.string for href in hrefs
            if href.string is not None and ".zip" in href.string]

    return pkgs


def list_packages(pkg_name=None):
    """
    List all packages, or all variants of a single package

    Parameters
    ----------
    pkg_name : str, optional
        - None: lists all stable packages on the server
        - <PackageName>: lists all variants of <PackageName> on the server

    Returns
    -------
    pkg_names : list

    Examples
    --------
    ::
        from scopesim import list_packages

        # list all stable packages on the server
        list_packages()

        # list all variants of a specific package
        list_packages("Armazones")

    """
    pkgs_dict = get_server_package_list()

    if pkg_name is None:
        pkg_names = list(pkgs_dict.keys())
    elif pkg_name in pkgs_dict:
        path = pkgs_dict[pkg_name]["path"]
        pkgs = get_server_folder_contents(path)
        pkg_names = [pkg for pkg in pkgs if pkg_name in pkg]

    return pkg_names


def download_packages(pkg_names, release="stable", save_dir=None, from_cache=None):
    """
    Download one or more packages to the local disk

    1. Download stable, dev
    2. Download specific version
    3. Download from github via url

    Parameters
    ----------
    pkg_names : str, list
        A list of package names, see ``list_packages()``

    release : str, optional
        By default, the most recent stable version of a package is downloaded.
        Other options are:
        - "stable" : the most recent stable version
        - "latest" : the latest development version to be published
        - a specific package filename as given by list_packages (see examples)
        - a github url for the specific branch and package (see examples)

    save_dir : str, optional
        The place on the local disk where the ``.zip`` package is to be saved.
        If left as None, defaults to the value in
        scopesim.rc.__config__["!SIM.file.local_packages_path"]

    from_cache : bool, optional
        Use the cached versions of the packages. If None, defaults to the RC
        value: ``!SIM.file.use_cached_downloads``

    Returns
    -------
    save_path : str
        The absolute path to the saved ``.zip`` package

    Examples
    --------
    ::
        from scopesim import download_packages, list_packages

        # Stable release of a list of packages
        download_packages(["test_package", "test_package"])

        # Development release of a single package
        download_packages("test_package", release="latest")

        # Specific version of the package
        list_packages("test_package")
        download_packages("test_package", release="2022-04-09.dev")

        # Specific package from a Gtihub commit hash or branch/tag name (use "@" or ":")
        download_packages("ELT", release="github:728761fc76adb548696205139e4e9a4260401dfc")
        download_packages("ELT", release="github@728761fc76adb548696205139e4e9a4260401dfc")
        download_packages("ELT", release="github@dev_master")

    """
    base_url = rc.__config__["!SIM.file.server_base_url"]

    pkgs_dict = get_server_package_list()

    if isinstance(pkg_names, str):
        pkg_names = [pkg_names]

    save_paths = []
    for pkg_name in pkg_names:
        if pkg_name in pkgs_dict:
            pkg_dict = pkgs_dict[pkg_name]
            path = pkg_dict["path"] + "/"

            from_github = False
            if release in ["stable", "latest"]:
                zip_name = pkg_dict[release]
                pkg_url = f"{base_url}{path}/{zip_name}.zip"
            elif "github" in release:
                base_url = "https://github.com/AstarVienna/irdb/tree/"
                github_hash = release.split(":")[-1].split("@")[-1]
                pkg_url = f"{base_url}{github_hash}/{pkg_name}"
                from_github = True
            else:
                zip_name = f"{pkg_name}.{release}.zip"
                pkg_variants = get_server_folder_contents(path)
                if zip_name not in pkg_variants:
                    raise ValueError(f"{zip_name} is not amoung the hosted "
                                     f"variants: {pkg_variants}")
                pkg_url = f"{base_url}{path}/{zip_name}"

            if save_dir is None:
                save_dir = rc.__config__["!SIM.file.local_packages_path"]
            if not os.path.exists(save_dir):
                os.mkdir(save_dir)

            if not from_github:
                try:
                    if from_cache is None:
                        from_cache = rc.__config__["!SIM.file.use_cached_downloads"]
                    cache_path = download_file(pkg_url, cache=from_cache)
                    save_path = os.path.join(save_dir, f"{pkg_name}.zip")
                    file_path = shutil.copy2(cache_path, save_path)

                    with zipfile.ZipFile(file_path, 'r') as zip_ref:
                        zip_ref.extractall(save_dir)

                except HTTPError as error:
                    raise ValueError(f"Unable to find file: {url + pkg_path}") from error
            else:
                download_github_folder(repo_url=pkg_url, output_dir=save_dir)
                save_path = save_dir

            save_paths += [os.path.abspath(save_path)]

        else:
            raise HTTPError(f"Unable to find package: {base_url + pkg_name}")

    return save_paths

# ==============================================================================
# Funtions below from from OLD_database.py
# ==============================================================================

# for backwards compatibility
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
    # todo: add proper depreciation warning
    text = "Function Depreciated --> please use scopesim.download_package-s-()"
    logging.warning(text)
    print(text)

    if isinstance(pkg_path, str):
        pkg_path = [pkg_path]

    pkg_names = [pkg.replace(".zip", "").split("/")[-1] for pkg in pkg_path]
    return download_packages(pkg_names, release="stable", save_dir=save_dir,
                             from_cache=from_cache)

def get_server_elements(url, unique_str="/"):
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


def list_example_data(url=None, return_files=False, silent=False):
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


def download_example_data(file_path, save_dir=None, url=None, from_cache=None):
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
    save_path : str
        The absolute path to the saved files
    """
    if isinstance(file_path, (list, tuple)):
        save_path = [download_example_data(thefile, save_dir, url)
                     for thefile in file_path]
    elif isinstance(file_path, str):

        if url is None:
            url = rc.__config__["!SIM.file.server_base_url"]
        if save_dir is None:
            save_dir = os.getcwd()
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        try:
            if from_cache is None:
                from_cache = rc.__config__["!SIM.file.use_cached_downloads"]
            cache_path = download_file(url + "example_data/" + file_path,
                                       cache=from_cache)
            save_path = os.path.join(save_dir, os.path.basename(file_path))
            file_path = shutil.copy2(cache_path, save_path)
        except HTTPError:
            ValueError(f"Unable to find file: {url + 'example_data/' + file_path}")

        save_path = os.path.abspath(save_path)

    return save_path


# """
# 2022-04-10 (KL)
# Code taken directly from https://github.com/sdushantha/gitdir
# Adapted for ScopeSim usage.
# Many thanks to the authors!
# """

def create_github_url(url):
    """
    From the given url, produce a URL that is compatible with Github's REST API. Can handle blob or tree paths.
    """
    repo_only_url = re.compile(r"https:\/\/github\.com\/[a-z\d](?:[a-z\d]|-(?=[a-z\d])){0,38}\/[a-zA-Z0-9]+$")
    re_branch = re.compile("/(tree|blob)/(.+?)/")

    # Check if the given url is a url to a GitHub repo. If it is, tell the
    # user to use 'git clone' to download it
    if re.match(repo_only_url,url):
        message = "✘ The given url is a complete repository. Use 'git clone' to download the repository"
        logging.error(message)
        raise ValueError(message)

    # extract the branch name from the given url (e.g master)
    branch = re_branch.search(url)
    download_dirs = url[branch.end():]
    api_url = (url[:branch.start()].replace("github.com", "api.github.com/repos", 1) +
              "/contents/" + download_dirs + "?ref=" + branch.group(2))
    return api_url, download_dirs


def download_github_folder(repo_url, output_dir="./"):
    """
    Downloads the files and directories in repo_url.

    Re-written based on the on the download function `here <https://github.com/sdushantha/gitdir/blob/f47ce9d85ee29f8612ce5ae804560a12b803ddf3/gitdir/gitdir.py#L55>`_
    """
    # convert repo_url into an api_url
    api_url, download_dirs = create_github_url(repo_url)

    # get the contents of the github folder
    user_interrupt_text = "GitHub download interrupted by User"
    try:
        opener = urllib.request.build_opener()
        opener.addheaders = [('User-agent', 'Mozilla/5.0')]
        urllib.request.install_opener(opener)
        response = urllib.request.urlretrieve(api_url)
    except KeyboardInterrupt:
        # when CTRL+C is pressed during the execution of this script
        logging.error(user_interrupt_text)
        raise ValueError(user_interrupt_text)

    # Make the base directories for this GitHub folder
    os.makedirs(os.path.join(output_dir, download_dirs), exist_ok=True)

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
                    opener.addheaders = [('User-agent', 'Mozilla/5.0')]
                    urllib.request.install_opener(opener)
                    # download the file
                    save_path = os.path.join(output_dir, entry['path'])
                    urllib.request.urlretrieve(entry["download_url"], save_path)
                    logging.info(f"Downloaded: {entry['path']}")

                except KeyboardInterrupt:
                    # when CTRL+C is pressed during the execution of this script
                    logging.error(user_interrupt_text)
                    raise ValueError(user_interrupt_text)
