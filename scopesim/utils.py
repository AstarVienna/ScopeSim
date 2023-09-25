"""Helper functions for ScopeSim."""
import math
from pathlib import Path
import sys
import logging
from collections import OrderedDict
from collections.abc import Iterable, Generator
from copy import deepcopy

from docutils.core import publish_string
import requests
import yaml
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.table import Column, Table

from . import rc


def unify(x, unit, length=1):
    """
    Convert all types of input to an astropy array/unit pair.

    Parameters
    ----------
    x : int, float, np.ndarray, astropy.Quantity
        The array to be turned into an ``astropy.Quantity``
    unit : astropy.Quantity
        The units to attach to the array
    length : int, optional
        If `x` is a scalar, and the desired output is an array with `length`.

    Returns
    -------
    y : astropy.Quantity
    """
    if isinstance(x, u.quantity.Quantity):
        if isinstance(x.value, np.ndarray):
            y = x.to(unit)
        elif length == 1:
            y = x.to(unit)
        else:
            y = ([x.value] * length * x.unit).to(unit)
    else:
        if hasattr(x, "__len__"):
            y = x * unit
        elif length == 1:
            y = x * unit
        else:
            y = [x] * length * unit

    return y


def parallactic_angle(ha, de, lat=-24.589167):
    r"""
    Compute the parallactic angle.

    Parameters
    ----------
    ha : float
        [hours] hour angle of target point
    de : float
        [deg] declination of target point
    lat : float
        [deg] latitude of observatory, defaults to Armazones

    Returns
    -------
    parang : float
       The parallactic angle

    Notes
    -----
    The parallactic angle is defined as the angle PTZ, where P is the
    .. math::
    \tan\eta = \frac{\cos\phi\sin H}{\sin\phi \cos\delta - \cos\phi \sin\delta \cos H}
    It is negative (positive) if target point is east (west) of the meridian.

    References
    ----------
    R. Ball: "A Treatise on Spherical Astronomy", Cambridge 1908
    """
    # Convert angles to radians
    ha = ha / 12. * np.pi
    de = np.deg2rad(de)
    lat = np.deg2rad(lat)

    eta = np.arctan2(np.cos(lat) * np.sin(ha),
                     np.sin(lat) * np.cos(de) -
                     np.cos(lat) * np.sin(de) * np.cos(ha))

    return np.rad2deg(eta)


def moffat(r, alpha, beta):
    """
    !!Unfinished!! Return a Moffat function.

    Parameters
    ----------
    r
    alpha
    beta

    Returns
    -------
    eta
    """
    return ((beta - 1) / (np.pi * alpha ** 2) *
            (1 + (r / alpha) ** 2) ** (-beta))


def poissonify(arr):
    """
    Add a realisation of the poisson process to the array `arr`.

    Parameters
    ----------
    arr : np.ndarray
        The input array which needs a Poisson distribution applied to items

    Returns
    -------
    arr : np.ndarray
        The input array, but with every pixel altered according to a poisson
        distribution
    """
    return np.random.poisson(arr).astype(np.float32)


def nearest(arr, val):
    """
    Return the index of the value from `arr` which is closest to `val`.

    Parameters
    ----------
    arr : np.ndarray, list, tuple
        Array to be searched
    val : float, int
        Value to find in `arr`

    Returns
    -------
    i : int
        index of array where the nearest value to `val` is
    """
    if isinstance(val, (list, tuple, np.ndarray)):
        arr = np.array(arr)
        return [nearest(arr, i) for i in val]

    return np.argmin(abs(arr - val))


def power_vector(val, degree):
    """Return the vector of powers of val up to a degree."""
    if degree < 0 or not isinstance(degree, int):
        raise ValueError("degree must be a positive integer")

    return np.array([val ** exp for exp in range(degree + 1)])


def deriv_polynomial2d(poly):
    """Derive (gradient) of a Polynomial2D model.

    Parameters
    ----------
    poly : astropy.modeling.models.Polynomial2D

    Returns
    -------
    gradient : tuple of Polynomial2d
    """
    import re
    from astropy.modeling.models import Polynomial2D
    degree = poly.degree
    dpoly_dx = Polynomial2D(degree=degree - 1)
    dpoly_dy = Polynomial2D(degree=degree - 1)
    regexp = re.compile(r"c(\d+)_(\d+)")
    for pname in poly.param_names:
        # analyse the name
        match = regexp.match(pname)
        i = int(match.group(1))
        j = int(match.group(2))
        cij = getattr(poly, pname)
        pname_x = "c%d_%d" % (i - 1, j)
        pname_y = "c%d_%d" % (i, j - 1)
        setattr(dpoly_dx, pname_x, i * cij)
        setattr(dpoly_dy, pname_y, j * cij)

    return dpoly_dx, dpoly_dy


def add_keyword(filename, keyword, value, comment="", ext=0):
    """
    Add a keyword, value pair to an extension header in a FITS file.

    Parameters
    ----------
    filename : str
        Name of the FITS file to add the keyword to
    keyword : str
    value : str, float, int
    comment : str
    ext : int, optional
        The fits extension index where the keyword should be added.
        Default is 0
    """
    f = fits.open(filename, mode="update")
    f[ext].header[keyword] = (value, comment)
    f.flush()
    f.close()


def airmass_to_zenith_dist(airmass):
    """
    Return zenith distance in degrees.

    Z = arccos(1/X)
    """
    return np.rad2deg(np.arccos(1. / airmass))


def zenith_dist_to_airmass(zenith_dist):
    """
    `zenith_dist` is in degrees.

    X = sec(Z)
    """
    return 1. / np.cos(np.deg2rad(zenith_dist))


def seq(start, stop, step=1):
    """Replacement for numpy.arange modelled after R's seq function.

    Returns an evenly spaced sequence from start to stop. stop is included if
    the difference between start and stop is an integer multiple of step.

    From the documentation of numpy.range: "When using a non-integer step, such
    as 0.1, the results will often not be consistent." This replacement aims to
    avoid these inconsistencies.

    Parameters
    ----------
    start, stop: [int, float]
        the starting and (maximal) end values of the sequence.

    step : [int, float]
        increment of the sequence, defaults to 1

    """
    feps = 1e-10  # value used in R seq.default

    delta = stop - start
    if delta == 0 and stop == 0:
        return stop
    try:
        npts = delta / step
    except ZeroDivisionError:
        if step == 0 and delta == 0:
            return start
        raise ValueError("invalid '(stop - start) / step'")

    if npts < 0:
        raise ValueError("wrong sign in 'step' argument")
    if npts > sys.maxsize:
        raise ValueError("'step' argument is much too small")

    reldd = abs(delta) / max(abs(stop), abs(start))

    if reldd < 100 * sys.float_info.epsilon:
        return start

    if isinstance(delta, int) and isinstance(step, int):
        # integer sequence
        npts = int(npts)
        return start + np.asarray(range(npts + 1)) * step

    npts = int(npts + feps)
    sequence = start + np.asarray(range(npts + 1)) * step
    # correct for possible overshot because of fuzz (from seq.R)
    if step > 0:
        return np.minimum(sequence, stop)
    else:
        return np.maximum(sequence, stop)


def add_mags(mags):
    """Return a combined magnitude for a group of py_objects with `mags`."""
    return -2.5 * np.log10((10 ** (-0.4 * np.array(mags))).sum())


def dist_mod_from_distance(d):
    """Use mu = 5 * np.log10(d) - 5 formula."""
    mu = 5 * np.log10(d) - 5
    return mu


def distance_from_dist_mod(mu):
    """Use d = 10**(1 + mu / 5) formula."""
    d = 10 ** (1 + mu / 5)
    return d


def telescope_diffraction_limit(aperture_size, wavelength, distance=None):
    """
    Return the diffraction limit of a telescope.

    Parameters
    ----------
    aperture_size : float
        [m] The diameter of the primary mirror

    wavelength : float
        [um] The wavelength for diffarction

    distance : float, optional
        Default is None. If ``distance`` is given, the transverse distance for
        the diffraction limit is returned in the same units as ``distance``


    Returns
    -------
    diff_limit : float
        [arcsec] The angular diffraction limit.
        If distance is not None, diff_limit is in the same units as distance

    """
    diff_limit = (((wavelength * u.um) / (aperture_size * u.m)) *
                  u.rad).to(u.arcsec).value

    if distance is not None:
        diff_limit *= distance / u.pc.to(u.AU)

    return diff_limit


def transverse_distance(angle, distance):
    """
    Turn an angular distance into a proper transverse distance.

    Parameters
    ----------
    angle : float
        [arcsec] The on-sky angle

    distance : float
        The distance to the object. Units are arbitary

    Returns
    -------
    trans_distance : float
        proper transverse distance. Has the same Units as ``distance``

    """
    trans_distance = angle * distance * u.AU.to(u.pc)

    return trans_distance


def angle_in_arcseconds(distance, width):
    """
    Return the angular distance of an object in arcseconds.

    Units must be consistent!
    """
    return np.arctan2(width, distance) * u.rad.to(u.arcsec)


def setup_loggers(**kwargs):
    """
    Set up both console and file loggers.

    Acceptable parameters are the same as the ``!SIM.logging`` sub dictionary

    """
    logd = rc.__currsys__["!SIM.logging"]
    logd.update(kwargs)

    logger = logging.getLogger()
    hdlr_names = [hdlr.name for hdlr in logger.handlers]

    if logd["log_to_file"] and "scopesim_file_logger" not in hdlr_names:
        f_handler = logging.FileHandler(logd["file_path"],
                                        logd["file_open_mode"])
        f_handler.name = "scopesim_file_logger"
        f_handler.setLevel(logd["file_level"])
        logger.addHandler(f_handler)

    if logd["log_to_console"] and "scopesim_console_logger" not in hdlr_names:
        s_handler = logging.StreamHandler(sys.stdout)
        s_handler.name = "scopesim_console_logger"
        s_handler.setLevel(logd["console_level"])
        logger.addHandler(s_handler)


def set_logger_level(which="console", level="ERROR"):
    """
    Set the level of logging for either the console or file logger.

    Parameters
    ----------
    which : {"console", "file"}
    level : {"ON", "OFF", "DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"}
    """
    hdlr_name = f"scopesim_{which}_logger"
    level = {"ON": "INFO", "OFF": "CRITICAL"}.get(level.upper(), level)
    logger = logging.getLogger()
    logger.setLevel(level)
    for hdlr in logger.handlers:
        if hdlr.name == hdlr_name:
            hdlr.setLevel(level)


def bug_report():
    """Get versions of dependencies for inclusion in bug report."""
    try:
        from importlib import import_module
    except ImportError:
        import_module = __import__

    packages = ["scopesim", "numpy", "scipy", "astropy", "matplotlib",
                "synphot", "skycalc_ipy", "requests", "bs4", "yaml"]

    # Check Python version
    print("Python:\n", sys.version)
    print("")

    # Check package dependencies
    for package_name in packages:
        try:
            pkg = import_module(package_name)
            print(package_name, ": ", pkg.__version__)
        except ImportError:
            print(package_name, "could not be loaded.")
        except AttributeError:
            print(package_name, ": version number not available")

    # Check operating system
    import platform
    osinfo = platform.uname()
    print("")
    print("Operating system: ", osinfo.system)
    print("         Release: ", osinfo.release)
    print("         Version: ", osinfo.version)
    print("         Machine: ", osinfo.machine)


def find_file(filename, path=None, silent=False):
    """Find a file in search path.

    Parameters
    ----------
    filename : str
        name of a file to look for
    path : list
        list of directories to search (default: ['./'])
    silent : bool
        if True, remain silent when file is not found

    Returns
    -------
    Absolute path of the file
    """
    if filename is None or filename.lower() == "none":
        return None

    if filename.startswith("!"):
        filename = from_currsys(filename)
    # Turn into pathlib.Path object for better manipulation afterwards
    filename = Path(filename)

    if path is None:
        path = rc.__search_path__

    if filename.is_absolute():
        # absolute path: only path to try
        trynames = [filename]
    else:
        # try to find the file in a search path
        trynames = [Path(trydir, filename)
                    for trydir in path if trydir is not None]

    for fname in trynames:
        if fname.exists():  # success
            # strip leading ./
            # Path should take care of this automatically!
            # while fname[:2] == './':
            #     fname = fname[2:]
            # Nevertheless, make sure this is actually the case...
            assert not str(fname).startswith("./")
            # HACK: Turn Path object back into string, because not everything
            #       that depends on this function can handle Path objects (yet)
            return str(fname)

    # no file found
    msg = f"File cannot be found: {filename}"
    if not silent:
        logging.error(msg)

    if from_currsys("!SIM.file.error_on_missing_file") is True:
        raise ValueError(msg)

    return None


def zendist2airmass(zendist):
    """Convert zenith distance to airmass.

    Parameters
    ----------
    zenith distance : [deg]
       Zenith distance angle

    Returns
    -------
    airmass in sec(z) approximation
    """
    return 1. / np.cos(np.deg2rad(zendist))


def airmass2zendist(airmass):
    """Convert airmass to zenith distance.

    Parameters
    ----------
    airmass : float (>= 1)

    Returns
    -------
    zenith distance in degrees
    """
    return np.rad2deg(np.arccos(1 / airmass))


def convert_table_comments_to_dict(tbl):

    comments_dict = {}
    if "comments" in tbl.meta:
        try:
            comments_str = "\n".join(tbl.meta["comments"])
            comments_dict = yaml.full_load(comments_str)
        except yaml.error.YAMLError:
            logging.warning("Couldn't convert <table>.meta['comments'] to dict")
            comments_dict = tbl.meta["comments"]
    elif "COMMENT" in tbl.meta:
        try:
            comments_dict = yaml.full_load("\n".join(tbl.meta["COMMENT"]))
        except yaml.error.YAMLError:
            logging.warning("Couldn't convert <table>.meta['COMMENT'] to dict")
            comments_dict = tbl.meta["COMMENT"]
    else:
        logging.debug("No comments in table")

    return comments_dict


def change_table_entry(tbl, col_name, new_val, old_val=None, position=None):

    offending_col = list(tbl[col_name].data)

    if old_val is not None:
        for ii in np.where(old_val in offending_col)[0]:
            offending_col[ii] = new_val
    elif position is not None:
        offending_col[position] = new_val
    else:
        raise ValueError("Either old_val or position must be given")

    fixed_col = Column(name=col_name, data=offending_col)

    ii = np.where(np.array(tbl.colnames) == col_name)[0][0]
    tbl.remove_column(col_name)
    tbl.add_column(fixed_col, index=ii)

    return tbl


def real_colname(name, colnames, silent=True):
    names = [name.lower(), name.upper(), name.capitalize()]
    real_name = [name for name in names if name in colnames]
    if not real_name:
        real_name = None
        if not silent:
            logging.warning("None of %s were found in %s", names, colnames)
    else:
        real_name = real_name[0]

    return real_name


def insert_into_ordereddict(dic, new_entry, pos):
    if isinstance(new_entry, dict):
        new_entry = [[key, val] for key, val in new_entry.items()]
    elif isinstance(new_entry, (list, tuple)) and \
            not isinstance(new_entry[0], (list, tuple)):
        new_entry = [new_entry]

    if pos < 0:
        pos += len(dic) + len(new_entry)

    new_dic = list(OrderedDict(dic).items())
    new_dic = new_dic[:pos] + new_entry + new_dic[pos:]
    new_dic = OrderedDict(new_dic)

    return new_dic


def empty_type(x):
    type_dict = {
        int: 0, float: 0., bool: False, str: " ",
        list: [], tuple: (), dict: {}
    }
    if "<U" in str(x):
        x = str

    return type_dict[x]


def get_meta_quantity(meta_dict, name, fallback_unit=""):
    """
    Extract a Quantity from a dictionary.

    Parameters
    ----------
    meta_dict : dict
    name : str
    fallback_unit : Quantity

    Returns
    -------
    quant : Quantity

    """
    if isinstance(meta_dict[name], str) and meta_dict[name].startswith("!"):
        meta_dict[name] = from_currsys(meta_dict[name])

    if isinstance(meta_dict[name], u.Quantity):
        unit = meta_dict[name].unit
    elif f"{name}_unit" in meta_dict:
        unit = meta_dict[f"{name}_unit"]
    else:
        unit = u.Unit(fallback_unit)

    quant = quantify(meta_dict[name], unit)

    return quant


def quantify(item, unit):
    """
    Ensure an item is a Quantity.

    Parameters
    ----------
    item : int, float, array, list, Quantity
    unit : str, Unit

    Returns
    -------
    quant : Quantity

    """
    if isinstance(item, str) and item.startswith("!"):
        item = from_currsys(item)
    if isinstance(item, u.Quantity):
        quant = item.to(u.Unit(unit))
    else:
        if isinstance(item, (np.ndarray, list, tuple)) and np.size(item) > 1000:
            quant = item << u.Unit(unit)
        else:
            quant = item * u.Unit(unit)
    return quant


def extract_type_from_unit(unit, unit_type):
    """
    Extract ``astropy`` physical type from a compound unit.

    Parameters
    ----------
    unit : astropy.Unit
    unit_type : str
        The physical type of the unit as given by ``astropy``

    Returns
    -------
    new_unit : Unit
        The input unit minus any base units corresponding to `unit_type`.
    extracted_units : Unit
        Any base units corresponding to `unit_type`.

    """
    extracted_units = u.Unit("")
    for base, power in zip(unit.bases, unit.powers):
        if unit_type == (base ** abs(power)).physical_type:
            extracted_units *= base ** power

    new_unit = unit / extracted_units

    return new_unit, extracted_units


def extract_base_from_unit(unit, base_unit):
    """
    Extract ``astropy`` base unit from a compound unit.

    Parameters
    ----------
    unit : astropy.Unit
    base_unit : Unit, str

   Returns
    -------
    new_unit : Unit
        The input unit minus any base units corresponding to `base_unit`.
    extracted_units : Unit
        Any base units corresponding to `base_unit`.

    """
    extracted_units = u.Unit("")
    for base, power in zip(unit.bases, unit.powers):
        if base == base_unit:
            extracted_units *= base ** power

    new_unit = unit * extracted_units ** -1

    return new_unit, extracted_units


def is_fits(filename) -> bool:
    # Using 'in ".fits"' to also catch ".fit", which exists sometimes...
    return (filename is not None and Path(filename).suffix.lower() in ".fits")


def get_fits_type(filename):
    with fits.open(filename) as hdulist:
        hdutype = "image"
        if hdulist[0].header["NAXIS"] == 0 and \
                hdulist[1].header["XTENSION"] == "BINTABLE":
            hdutype = "bintable"

    return hdutype


def quantity_from_table(colname: str, table: Table,
                        default_unit: str = "") -> u.Quantity:
    col = table[colname]
    if col.unit is not None:
        return col.quantity

    unit = unit_from_table(colname, table, default_unit)
    # TODO: or rather << ?
    return col * unit


def unit_from_table(colname: str, table: Table,
                    default_unit: str = "") -> u.Unit:
    """
    Look for the unit for a column based on the meta dict keyword "<col>_unit".
    """
    col = table[colname]
    if col.unit is not None:
        return col.unit

    colname_u = f"{colname}_unit"
    if colname_u in table.meta:
        return u.Unit(table.meta[colname_u])

    com_tbl = convert_table_comments_to_dict(table)
    if colname_u in com_tbl:
        return u.Unit(com_tbl[colname_u])

    tbl_name = table.meta.get("name", table.meta.get("filename"))
    logging.info(("%s_unit was not found in table.meta: %s. Default to: %s"),
                 colname, tbl_name, default_unit)

    return u.Unit(default_unit)


def deg2rad(theta):
    return theta * math.pi / 180


def rad2deg(theta):
    return theta * 180 / math.pi


def has_needed_keywords(header, suffix=""):
    """Check to see if the WCS keywords are in the header."""
    keys = {"CDELT1", "CRVAL1", "CRPIX1"}
    keys = {key + suffix for key in keys}
    keys.add("NAXIS1")
    return all(key in header.keys() for key in keys)


def stringify_dict(dic, ignore_types=(str, int, float)):
    """Turn a dict entries into strings for addition to FITS headers."""
    dic_new = deepcopy(dic)
    for key in dic_new:
        if not isinstance(dic_new[key], ignore_types):
            dic_new[key] = str(dic_new[key])

    return dic_new


def clean_dict(orig_dict, new_entries):
    """
    Replace OBS_DICT keywords with actual values.

    Parameters
    ----------
    orig_dict : dict

    new_entries : dict
        OBS dict

    Returns
    -------
    orig_dict : dict
        Updated dict

    """
    for key in orig_dict:
        if isinstance(orig_dict[key], str) and orig_dict[key] in new_entries:
            orig_dict[key] = new_entries[orig_dict[key]]

    return orig_dict


def from_currsys(item):
    """Return the current value of a bang-string from ``rc.__currsys__``."""
    if isinstance(item, Table):
        tbl_dict = {col: item[col].data for col in item.colnames}
        tbl_dict = from_currsys(tbl_dict)
        item_meta = item.meta
        item = Table(data=list(tbl_dict.values()),
                     names=list(tbl_dict.keys()))
        item.meta = item_meta

    if isinstance(item, np.ndarray) and not isinstance(item, u.Quantity):
        item = np.array([from_currsys(x) for x in item])

    if isinstance(item, list):
        item = [from_currsys(x) for x in item]

    if isinstance(item, dict):
        for key in item:
            item[key] = from_currsys(item[key])

    if isinstance(item, str) and len(item) and item.startswith("!"):
        if item in rc.__currsys__:
            item = rc.__currsys__[item]
        else:
            raise ValueError(f"{item} was not found in rc.__currsys__")

    if isinstance(item, str) and item.lower() == "none":
        item = None

    return item


def check_keys(input_dict, required_keys, action="error", all_any="all"):
    """Check to see if all/any of the required keys are present in a dict."""
    if isinstance(input_dict, (list, tuple)):
        input_dict = {key: None for key in input_dict}

    if all_any == "all":
        keys_present = all(key in input_dict for key in required_keys)
    elif all_any == "any":
        keys_present = any(key in input_dict for key in required_keys)
    else:
        raise ValueError("all_any must be either 'all' or 'any'")

    if not keys_present:
        if "error" in action:
            raise ValueError("One or more of the following keys missing from "
                             f"input_dict: \n{required_keys} "
                             f"\n{input_dict.keys()}")
        if "warn" in action:
            logging.warning(("One or more of the following keys missing "
                             "from input_dict: \n%s \n%s"), required_keys,
                            input_dict.keys())

    return keys_present


def interp2(x_new, x_orig, y_orig):
    """Check and correct for decreasing x_orig values."""
    if x_orig[0] < x_orig[-1]:
        y_new = np.interp(x_new, x_orig, y_orig)
    else:
        y_new = np.interp(x_new, x_orig[::-1], y_orig[::-1])

    return y_new


def write_report(text, filename=None, output=None):
    """Write a report string to file in latex or rst format."""
    if output is None:
        output = ["rst"]
    elif isinstance(output, str):
        output = [output]

    if filename is not None:
        for fmt in output:
            out_text = deepcopy(text)
            if fmt.lower() == "latex":
                out_text = publish_string(out_text, writer_name="latex")
                out_text = out_text.decode("utf-8")

            suffix = {"rst": ".rst", "latex": ".tex"}[fmt]
            fname = Path(filename).with_suffix(suffix)
            fname.write_text(out_text, encoding="utf-8")


def pretty_print_dict(dic, indent=0):
    text = ""
    for key, value in dic.items():
        if isinstance(value, dict):
            text += " " * indent + f"{str(key)}:\n"
            text += pretty_print_dict(value, indent=indent + 2)
        else:
            text += " " * indent + f"{str(key)}: {str(value)}\n"

    return text


def return_latest_github_actions_jobs_status(
        owner_name="AstarVienna",
        repo_name="ScopeSim",
        branch="dev_master",
        actions_yaml_name="tests.yml",
    ):
    """Get the status of the latest test run."""
    response = requests.get(
        f"https://api.github.com/repos/{owner_name}/{repo_name}/actions/"
        f"workflows/{actions_yaml_name}/runs?branch={branch}&per_page=1"
    )
    dic = response.json()
    run_id = dic["workflow_runs"][0]["id"]

    response = requests.get(
        f"https://api.github.com/repos/{owner_name}/{repo_name}/actions/runs/"
        f"{run_id}/jobs"
    )
    dic = response.json()
    params_list = []
    for job in dic["jobs"]:
        params = {
            "name": job["name"],
            "status": job["status"],
            "conclusion": job["conclusion"],
            "started_at": job["started_at"],
            "completed_at": job["completed_at"],
            "url": job["html_url"],
            "badge_url": None
        }

        # TODO: this could use the new badges from IRDB, once that's in
        #       scopesim_core...
        key = "Python_" + job["name"].split()[-1][:-1]
        value = "passing" if job["conclusion"] == "success" else "failing"
        colour = "brightgreen" if job["conclusion"] == "success" else "red"
        badge_url = f"https://img.shields.io/badge/{key}-{value}-{colour}"
        params["badge_url"] = badge_url
        params_list.append(params)

    return params_list


def close_loop(iterable: Iterable) -> Generator:
    """x, y = zip(*close_loop(zip(x, y)))"""
    iterator = iter(iterable)
    first = next(iterator)
    yield first
    yield from iterator
    yield first


def figure_factory(nrows=1, ncols=1, **kwargs):
    """Default way to init fig and ax, to easily modify later."""
    iterable_axes = kwargs.pop("iterable_axes", False)
    fig, ax = plt.subplots(nrows, ncols, **kwargs)
    if iterable_axes and not isinstance(ax, Iterable):
        ax = (ax,)
    return fig, ax


def figure_grid_factory(nrows=1, ncols=1, **kwargs):
    """Gridspec variant."""
    fig = plt.figure()
    gs = fig.add_gridspec(nrows, ncols, **kwargs)
    return fig, gs
