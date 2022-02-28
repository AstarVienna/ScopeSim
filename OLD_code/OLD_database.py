"""
Put functions which connect the IPDB to SimCADO here
Possibly also add the skycalc_cli interface here

Write functions to:
0. Connect to the server
1a. Read which packages are available
1b. Look at which packages are available locally
2. Display a list of which packages are available
3. Download a package
4. Unpack into it's own folder

"""

import os
import requests
import datetime as dt
import zipfile as zf
import logging

from numpy import where as npwhere
from astropy.table import Table, Row, vstack
from astropy.io import ascii as ioascii

from scopesim.utils import download_file
from scopesim import rc

__all__ = ["list_all", "list_instruments", "list_psfs", "list_source_packages",
           "get_local_packages", "get_server_packages",
           "download_package", "set_up_local_package_directory",
           "local_db_paths", "server_db_urls",
           "find_package_on_disk", "find_package_on_server", "get_path"]


LOCAL_DB_HEADER_PATTERN = """# Date-created : {}
# Date-modified : {}
# Description : Packages containing {} specific data files
#
name   author   date_added  date_modified   path     
"""


def _local_inst_db():
    return os.path.join(rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                        rc.__rc__["FILE_INST_PKG_LOCAL_DB_NAME"])


def _local_psf_db():
    return os.path.join(rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                        rc.__rc__["FILE_PSF_LOCAL_DB_NAME"])


def _local_src_db():
    return os.path.join(rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                        rc.__rc__["FILE_SRC_PKG_LOCAL_DB_NAME"])


def set_local_path_names(path):
    local_inst_db = os.path.join(path, rc.__rc__["FILE_INST_PKG_LOCAL_DB_NAME"])
    local_psf_db = os.path.join(path, rc.__rc__["FILE_PSF_LOCAL_DB_NAME"])
    local_src_db = os.path.join(path, rc.__rc__["FILE_SRC_PKG_LOCAL_DB_NAME"])

    return local_inst_db, local_psf_db, local_src_db


def _local_paths():
    return set_local_path_names(rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"])


def _svr_inst_db():
    return rc.__rc__["FILE_SERVER_BASE_URL"] + \
           rc.__rc__["FILE_INST_PKG_SERVER_DB_NAME"]


def _svr_psf_db():
    return rc.__rc__["FILE_SERVER_BASE_URL"] + \
           rc.__rc__["FILE_PSF_SERVER_DB_NAME"]


def _svr_src_db():
    return rc.__rc__["FILE_SERVER_BASE_URL"] + \
           rc.__rc__["FILE_SRC_PKG_SERVER_DB_NAME"]


def local_db_paths(name=None):
    """
    Return the paths for the local database files

    Parameters
    ----------
    name : str, optional
        None, "inst", "psf", "st"

    Returns
    -------
    svr_dict : dict or str

    """

    local_dict = {"inst" : _local_inst_db(), "psf" : _local_psf_db(),
                  "st" : _local_src_db()}
    if name is None:
        return local_dict
    else:
        return local_dict[name]


def server_db_urls(name=None):
    """
    Return the URLs for the server side database files

    Parameters
    ----------
    name : str, optional
        None, "inst", "psf", "st"

    Returns
    -------
    svr_dict : dict or str

    """

    svr_dict = {"inst": _svr_inst_db(), "psf": _svr_psf_db(),
                "st": _svr_src_db()}
    if name is None:
        return svr_dict
    else:
        return svr_dict[name]


def set_up_local_package_directory(dirname=None, overwrite=False):
    """
    Sets up the directory structure for a set of locally stored package files

    Parameters
    ----------
    dirname : str, optional
        Path to where you would like to store downloaded packages.
        If ``None``, the value from the SimCADO RC file is taken

    overwrite : bool, optional
        Default False. Whether to overwrite an existing directory structure


    Examples
    --------
    ::

        >>> import scopesim.server as svr
        >>> svr.set_up_local_package_directory("./scopesim_downloads/")


    See Also
    --------
    Using the SimCADO RC file

    """

    if dirname is None:
        dirname = rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"]

    # set up downloads directory and sub directories
    for dname in [dirname,
                  os.path.join(dirname, rc.__rc__["FILE_SCOPE_PKG_LOCAL_PATH"]),
                  os.path.join(dirname, rc.__rc__["FILE_INST_PKG_LOCAL_PATH"]),
                  os.path.join(dirname, rc.__rc__["FILE_PSF_LOCAL_PATH"]),
                  os.path.join(dirname, rc.__rc__["FILE_SRC_PKG_LOCAL_PATH"])]:
        if not os.path.exists(dname):
            os.makedirs(dname)
        elif not overwrite:
            print("{} already exists. If you would like to overwrite it, set"
                  "overwrite=True".format(dname))

    # make db files
    now = dt.datetime.now().strftime('%Y-%m-%d')
    local_paths = set_local_path_names(dirname)
    package_type = ["instrument/telescope", "PSF", "Source object"]

    for loc_db, pkg_type in zip(local_paths, package_type):
        if not os.path.exists(loc_db) or overwrite:
            with open(loc_db, "w") as new_db:
                hdr_text = LOCAL_DB_HEADER_PATTERN.format(now, now, pkg_type)
                new_db.write(hdr_text)
        elif not overwrite:
            print("{} already exists. If you would like to overwrite it, set"
                  "'overwrite=True'".format(loc_db))


def get_local_packages(path=None):
    """
    Return a list of packages on the local disk

    By default returns only the instrument package server database
    To get the PSF or source databases, pass relevant file path.
    These can be found under::

        scopesim.server.local_db_paths(<name>)

    where ``name`` can be ``None``,``inst``, ``psf``, or ``st``


    Parameters
    ----------
    path : str
        URL to package database file on the server. See above.

    Returns
    -------
    local_table : `astropy.Table`

    Examples
    --------
    ::

        >>> import scopesim.server as svr
        >>> psf_local_path = svr.local_db_paths("psf")
        >>> svr_pkgs = svr.get_server_packages(psf_local_path)

    """

    if path is None:
        path = _local_inst_db()

    if not os.path.exists(path):
        raise ValueError(path + " doesn't exist")

    # If this throws an error, it's because the DB file is not formatted
    # correctly
    # The column row names should NOT be commented out
    local_table = ioascii.read(path, format="basic", fast_reader=False)
    local_table = rename_table_colnames(local_table)

    # make sure all columns are strings
    for col in local_table.colnames:
        local_table[col] = local_table[col].astype("str")

    return local_table


def get_server_text(path=None):
    """Returns the text response from the server for a DB file"""

    if path is None:
        path = _svr_inst_db()

    try:
        server_db_text = requests.get(path).text
        server_db_text = server_db_text.replace("\r", "")
    except:
        logging.warning("Connection could not be established to "
                      "{}".format(_svr_inst_db()))
        server_db_text = None

    return server_db_text


def rename_table_colnames(svr_table):
    """Reinstated column names if they where accidentally parsed to comments"""
    if svr_table.colnames[0] != "col0":
        return svr_table

    col_names = svr_table.meta["comments"][-1].replace("#", "").split()

    if len(col_names) != len(svr_table.colnames):
        raise ValueError("Somethings up with the server database header names")

    for i, col_name in enumerate(col_names):
        svr_table[svr_table.colnames[i]].name = col_name

    return svr_table


def get_server_packages(path=None):
    """
    Return a list of packages on the server

    By default returns only the instrument package server database
    To get the PSF or source databases, pass relevant url. These can be found
    under::

        scopesim.server.svr_db_dict(<name>)

    where ``name`` can be ``None``,``inst``, ``psf``, or ``st``


    Parameters
    ----------
    path : str
        URL to package database file on the server. See above.

    Returns
    -------
    local_table : `astropy.Table`

    Examples
    --------
    ::

        >>> import scopesim.server as svr
        >>> psf_server_url = svr.server_db_urls("psf")
        >>> svr_pkgs = svr.get_server_packages(psf_server_url)

    """

    svr_text = get_server_text(path)
    if svr_text is None:
        return None
    svr_table = ioascii.read(svr_text)
    svr_table = rename_table_colnames(svr_table)

    return svr_table


def list_packages(local_path=None, server_url=None, return_table=False,
                  msg="Packages"):
    """
    Prints to screen (returns) a list of all available packages
        
    Parameters
    ----------
    local_path : str, optional
        Default `scopesim.server._local_psf_db()`

    server_url : str, optional
        Default `scopesim.server._svr_psf_db()`

    return_table : bool, optional
    msg : str

    Returns
    -------
    all_table : `astropy.Table`
        Only if `return_table` is `True`

    """

    local_table = get_local_packages(local_path)
    svr_table = get_server_packages(server_url)
    if svr_table is None:
        all_table = local_table
    else:
        all_table = vstack(local_table, svr_table)

    print("\n{} saved offline\n".format(msg) + "="*(len(msg)+14))
    print(local_table)
    print("\n{} on the server\n".format(msg) + "="*(len(msg)+14))
    print(svr_table)

    if return_table:
        return all_table


def list_instruments(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available instrument packages

    By default `list_instruments` looks in

    * `scopesim.server.local_db_paths("inst")`
    * `scopesim.server.server_db_urls("inst")`

    Parameters
    ----------
    local_path : str, optional
        Default `scopesim.server._local_psf_db()`

    server_url : str, optional
        Default `scopesim.server._svr_psf_db()`

    return_table : bool, optional
    msg : str

    Returns
    -------
    all_table : `astropy.Table`
        Only if `return_table` is `True`


    See Also
    --------
    :func:`.list_packages`

    """

    if local_path is None:
        local_path = _local_inst_db()
    if server_url is None:
        server_url = _svr_inst_db()

    return list_packages(local_path, server_url, return_table,
                         "Instrument and Telescope packages")


def list_psfs(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available PSF files

    By default `list_psfs` looks in

    * `scopesim.server.local_db_paths("psf")`
    * `scopesim.server.server_db_urls("psf")`

    Parameters
    ----------
    local_path : str, optional
        Default `scopesim.server._local_psf_db()`

    server_url : str, optional
        Default `scopesim.server._svr_psf_db()`

    return_table : bool, optional
    msg : str

    Returns
    -------
    all_table : `astropy.Table`
        Only if `return_table` is `True`


    See Also
    --------
    :func:`.list_packages`

    """

    if local_path is None:
        local_path = _local_psf_db()
    if server_url is None:
        server_url = _svr_psf_db()

    return list_packages(local_path, server_url, return_table,
                         "PSF files")


def list_source_packages(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available source packages

    By default `list_source_pkgs` looks in

    * `scopesim.server.local_db_paths("st")`
    * `scopesim.server.server_db_urls("st")`


    Parameters
    ----------
    local_path : str, optional
        Default `scopesim.server._local_psf_db()`

    server_url : str, optional
        Default `scopesim.server._svr_psf_db()`

    return_table : bool, optional
    msg : str


    Returns
    -------
    all_table : `astropy.Table`
        Only if `return_table` is `True`


    See Also
    --------
    :func:`.list_packages`

    """

    if local_path is None:
        local_path = _local_src_db()
    if server_url is None:
        server_url = _svr_src_db()

    return list_packages(local_path, server_url, return_table,
                         "Source packages")


def list_all():
    """
    Prints all the packages on the server and on the local disk to screen

    See Also
    --------
    :func:`.list_instruments`
    :func:`.list_psfs`
    :func:`.list_source_packages`

    """
    list_instruments()
    list_psfs()
    list_source_packages()


def check_package_exists(pkg_name, svr_path=None):
    """Checks on the server that a certain package name exists"""
    if svr_path is None:
        svr_path = _svr_inst_db()

    svr_base_url = os.path.dirname(svr_path)
    svr_table = get_server_packages(svr_path)
    pkg_entry = get_package_table_entry(pkg_name, svr_table)

    if pkg_entry is not None:
        url = os.path.join(svr_base_url, pkg_entry["path"]).replace("\\", "/")

        if not requests.get(url).ok:
            raise ValueError(url + " doesn't return ALL_GOOD (200)")

        return_val = True
    else:
        return_val = False

    return return_val


def get_server_package_path(pkg_name, svr_table):
    """Find the url to a package ZIP file on the server"""

    pkg_path = None
    if pkg_name in svr_table["name"]:
        svr_dict = {n: p for n, p in zip(svr_table["name"], svr_table["path"])}
        pkg_path = svr_dict[pkg_name]

    return pkg_path


def get_package_table_entry(pkg_name, db_table):
    """Returns an astropy.table.Row object with the package data"""

    # ::todo implement multiple entry handling
    if pkg_name in db_table["name"]:
        pkg_index = npwhere(db_table["name"] == pkg_name)[0]
        pkg_entry = db_table[pkg_index[0]]
    else:
        pkg_entry = None

    return pkg_entry


def determine_type_of_package(svr_db_filename):
    """Guess what the package type is: 'inst', 'psf', 'st'"""

    pkg_type = None
    if "inst" in svr_db_filename.lower():
        pkg_type = "inst"
    elif "psf" in svr_db_filename.lower():
        pkg_type = "psf"
    elif "source" in svr_db_filename.lower():
        pkg_type = "st"

    return pkg_type


def extract_package(pkg_name, overwrite=True):
    """
    Extracts the files from the ZIP archive for a downloaded instrument package

    Parameters
    ----------
    pkg_name : str
        the name of the zip file containing instrument package data

    overwrite : bool
        Flag to overwrite already extracted data. Only set to True if you are
        sure you don't mind losing any changes to the package data

    """

    if os.path.exists(pkg_name) and ".zip" in pkg_name:
        file_path = pkg_name
    else:
        pkg_entry = find_package_on_disk(pkg_name)
        if not isinstance(pkg_entry, Row):
            raise ValueError("{} wasn't found on disk".format(pkg_name))

        file_path = os.path.join(rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"],
                                 pkg_entry["path"])

    new_dir = file_path.replace(".zip", "")
    if os.path.exists(new_dir) and not overwrite:
        pass
    else:
        if os.path.exists(new_dir):
            print("{} exists, but overwriting anyway".format(pkg_name))
        with zf.ZipFile(file_path) as pkg_zip:
            pkg_zip.extractall(new_dir)


def download_package(pkg_name, unzip_package=True, save_dir=None,
                     server_dbs=None):
    """
    Download a package from the server and extract into the downloads folder

    .. warning::
        By default ``unzip_package=True``. This will overwrite an
        existing package with the same **FILENAME**, not of the same package
        name. (Files are named ``INSTRUMENT_YYYY-MM-DD.zip``). If you have
        altered any individual files in the package directory, these will be
        overwritten if you re-download the same package zip file.


    Parameters
    ----------
    pkg_name : str
        Name of the package to download

    unzip_package : bool, optional
        Default True. If True, the package is automatically extracted into a
        folder with the same name as the zip file. If False, only the zip file
        is downloaded

    save_dir : str, optional
        Where to save the package on the local disk. Default INST_PKG_LOCAL_PATH

    server_dbs : str, optional
        URL to the server


    Returns
    -------
    local_filename : str
        The path and filename of the saved zip file


    Examples
    --------
    ``download_package`` queries all the database files on the server and
    downloads the first package it finds which matched ``pkg_name``::

        >>> # Get the 'MICADO' package from the list of instruments
        >>> download_package("MICADO")
        >>> # Get the 'NIR_Sources' package from the list of source spectra
        >>> download_package("NIR_Sources")

    By default ``download_package`` extracts the zip file after downloading.
    If we just want to download the package and save it somewhere else, we can
    set the ``unzip_package`` flag to False::

        >>> download_package("MICADO", unzip_package=False,
        ...                  save_dir="./random/folder/")

    """

    pkg_entry, svr_db = find_package_on_server(pkg_name,
                                               server_dbs=server_dbs,
                                               return_db_filename=True)

    if pkg_entry is None:
        raise ValueError("{} wasn't found on the server".format(pkg_name))

    pkg_url = rc.__rc__["FILE_SERVER_BASE_URL"] + pkg_entry["path"]
    pkg_type = determine_type_of_package(svr_db)

    if not check_package_exists(pkg_name, server_db_urls()[pkg_type]):
        raise ValueError("Package is missing: " + pkg_name)

    if save_dir is None:
        stem = os.path.dirname(pkg_entry["path"])
        save_dir = os.path.join(rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"], stem)

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    local_filename = download_file(pkg_url, save_dir)
    print("Saved {} in {}".format(pkg_name, local_filename))
    if unzip_package and ".fits" not in local_filename:
        extract_package(local_filename, overwrite=True)
        unzip_dir = local_filename.replace(".zip", "")
        print("Unzipped {} to {}".format(pkg_name, unzip_dir))

    local_db_path = local_db_paths()[pkg_type]
    new_local_tbl = add_pkg_to_local_db(pkg_entry, local_db_path)

    write_table_to_disk(new_local_tbl, local_db_path)

    return local_filename


def write_table_to_disk(tbl, path):
    """Write an database file back to disk. Overwrite existing DB file"""
    tbl.write(path, format="ascii.fixed_width", overwrite=True, delimiter="")


def find_package_on_server(pkg_name, server_dbs=None, return_db_filename=False):
    """
    Returns the first match for ``pkg_name`` in the server database files

    Parameters
    ----------
    pkg_name : str
    server_dbs : list
    return_db_filename : bool

    Returns
    -------
    pkg_entry : :class:`astropy.table.Row`
        The row data for the package from the server database tabl

    svr_db_basename : str
        If ``return_db_filename==True``, the database file name, so that the
        type of package can be determined

    """

    if server_dbs is None:
        server_dbs = [_svr_inst_db(), _svr_psf_db(), _svr_src_db()]

    pkg_entry, svr_db = None, None
    for svr_db in server_dbs:
        svr_table = get_server_packages(svr_db)
        if svr_table is None:
            break
        pkg_entry = get_package_table_entry(pkg_name, svr_table)
        if pkg_entry is not None:
            break

    if return_db_filename:
        return pkg_entry, os.path.basename(svr_db)
    else:
        return pkg_entry


def find_package_on_disk(pkg_name, local_dbs=None, return_db_filename=False):
    """
    Returns the first match for ``pkg_name`` in the local database files

    Parameters
    ----------
    pkg_name : str
    local_dbs : list
    return_db_filename : bool

    Returns
    -------
    pkg_entry : :class:`astropy.Table`, None
        Returns ``None`` if package not found

    """

    if local_dbs is None:
        local_dbs = [_local_inst_db(), _local_psf_db(), _local_src_db()]

    pkg_entry, local_db = None, None
    for local_db in local_dbs:
        local_table = get_local_packages(local_db)
        pkg_entry = get_package_table_entry(pkg_name, local_table)
        if pkg_entry is not None:
            break

    if return_db_filename:
        return pkg_entry, os.path.basename(local_db)
    else:
        return pkg_entry


def add_pkg_to_local_db(new_row, local_db):
    """
    Adds a package entry to the local database file

    Parameters
    ----------
    new_row : ``astropy.table.Row``
        The package entry data

    local_db : str, ``astropy.Table``
        Eith


    Returns
    -------

    """

    if isinstance(local_db, str):
        local_table = get_local_packages(local_db)
    elif isinstance(local_db, Table):
        local_table = local_db
    else:
        raise ValueError("local_db must be either Table or path to DB file")

    if type(new_row) != Row:
        raise ValueError("pkg_entry must be an astropy.table.Row object")

    if new_row["name"] in local_table["name"]:
        ii = npwhere(local_table["name"] == new_row["name"])[0][0]
        fmt = '%Y-%m-%d'
        local_date = dt.datetime.strptime(local_table[ii]["date_modified"], fmt)
        new_date   = dt.datetime.strptime(new_row["date_modified"], fmt)

        if new_date >= local_date:

            dic = {col : local_table[ii][col] for col in local_table.colnames}
            dic["name"] = dic["name"] + "_" + dic["date_modified"]
            new_tbl = Table(names=[col for col in dic],
                            data=[[dic[col]] for col in dic])
            new_tbl_2 = Table(names=[col for col in new_row.colnames],
                              data=[[new_row[col]] for col in new_row.colnames])
            tbl_to_add = vstack([new_tbl, new_tbl_2])

            local_table.remove_row(ii)
        else:
            dic = {col: new_row[col] for col in new_row.colnames}
            dic["name"] = dic["name"] + "_" + dic["date_modified"]
            tbl_to_add = Table(names=[col for col in dic],
                               data=[[dic[col]] for col in dic])

    else:
        tbl_to_add = Table(names=[col for col in new_row.colnames],
                           data=[[new_row[col]] for col in new_row.colnames])

    new_local_table = vstack([local_table, tbl_to_add])
    new_local_table.meta = local_table.meta

    return new_local_table


def get_path(pkg_name):
    """Returns the absolute path of package ``pkg_name``"""

    abs_path = None
    pkg_entry = find_package_on_disk(pkg_name)
    if pkg_entry is not None:
        downloads_path = rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"]
        abs_path = os.path.join(downloads_path, pkg_entry["path"])

    return abs_path
