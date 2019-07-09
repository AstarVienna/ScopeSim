import requests
import bs4

from scopesim import rc


def get_server_elements(url=None, unique_str="/"):
    dirs = []
    if url is not None:
        result = requests.get(url).content
        soup = bs4.BeautifulSoup(result)
        dirs = soup.findAll("a", href=True)
        dirs = [tmp.string for tmp in dirs if unique_str in tmp.string]

    return dirs


def get_packages(url=None):
    if url is None:
        url = rc.__rc__["FILE_SERVER_BASE_URL"]

    all_pkgs = []
    folders = get_server_elements(url, "/")
    for folder in folders:
        pkgs = get_server_elements(url + folder, ".zip")
        all_pkgs += [folder + pkg for pkg in pkgs]

    return all_pkgs
