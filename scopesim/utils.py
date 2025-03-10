# -*- coding: utf-8 -*-
"""Helper functions for ScopeSim."""

from pathlib import Path
import sys
import logging
from logging.config import dictConfig
from collections.abc import Iterable, Generator, Set, Mapping
from copy import deepcopy
from typing import TextIO, Union
from io import StringIO
from importlib import metadata
import functools

from docutils.core import publish_string
import yaml
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.table import Column, Table

from astar_utils import get_logger, is_bangkey

from . import rc


logger = get_logger(__name__)
bug_logger = get_logger("bug_report")


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


# TODO: I think we have multiple implementations of such a thing across out
#       various packages. Should be harmonized and go into astar-utils.
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


def _get_required_packages():
    reqs = metadata.requires(__package__)
    # metadata.requires can return None if the package metadata cannot be found
    if reqs is None:
        return []
    for req in reqs:
        # Only include non-extra packages
        if "extra" in req:
            continue

        name = req.split(">", maxsplit=1)[0].strip()
        name = name.split("(", maxsplit=1)[0].strip()
        yield name


def _get_all_irdb_pkgs(root: Path):
    return [pkg_path for pkg_path in root.iterdir() if pkg_path.is_dir()
            and not pkg_path.name.startswith("__")] if root.is_dir() else []


def _get_irdb_pkg_version(pkg_path: Path) -> str:
    versionfile = pkg_path / "version.yaml"
    if not versionfile.exists():
        return "version number not available."
    with versionfile.open(encoding="utf-8") as file:
        return yaml.load(file, yaml.SafeLoader)["version"]


def _write_bug_report(stream: TextIO) -> None:
    # Check Python version
    stream.write(f"Python:\n{sys.version}\n")

    # Check package dependencies
    stream.write("\nInstalled Python packages:\n")
    packages = set(_get_required_packages())
    packages.update({"scopesim", "scopesim_templates", "scopesim_data", "anisocado"})
    maxkeylen = max(len(pkg) for pkg in packages)
    for package_name in sorted(packages):
        stream.write(f"{package_name:>{maxkeylen+2}}: ")
        try:
            ver = metadata.version(package_name)
            stream.write(f"{ver}\n")
        except ImportError:
            stream.write("could not be loaded.\n")
        # except AttributeError:
        #     stream.write(f"version number not available.\n")

    # Check IRDB packages
    stream.write("\nInstalled IRDB packages:\n")
    pkgs_path = Path(rc.__config__["!SIM.file.local_packages_path"])
    installed_pkgs = _get_all_irdb_pkgs(pkgs_path)
    maxkeylen = max((len(pkg.stem) for pkg in installed_pkgs), default=0)
    for pkg_path in installed_pkgs:
        pkg_ver = _get_irdb_pkg_version(pkg_path)
        stream.write(f"{pkg_path.stem:>{maxkeylen+2}}: {pkg_ver}\n")

    # Check operating system
    import platform
    osinfo = platform.uname()
    stream.write("\nOperating System info:\n")
    for field in ["system", "release", "version", "machine"]:
        stream.write(f"{field.title():>9}: {getattr(osinfo, field)}\n")


def bug_report() -> None:
    """Print versions of dependencies for inclusion in bug report."""
    _write_bug_report(sys.stdout)


def bug_report_to_file(filename) -> None:
    """Like bug_report, but writes to file instead of printing."""
    filename = Path(filename)
    with filename.open("w", encoding="utf-8") as file:
        _write_bug_report(file)


def log_bug_report(level=logging.DEBUG) -> None:
    """Emit bug report as logging message."""
    with StringIO() as str_stream:
        _write_bug_report(str_stream)
        bug_logger.log(level, str_stream.getvalue())


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
        raise ValueError(f"!-string filename should be resolved upstream: "
                         f"{filename}")
        # filename = from_currsys(filename)
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
        logger.error(msg)

    # TODO: Not sure what to do here
    if from_currsys("!SIM.file.error_on_missing_file"):
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


# TODO: There are identical implementations of these functions with slightly
#       different names in this module. The ones WITHOUT underscores are
#       actually used, better documented and have unit tests. But the one WITH
#       underscores have more readable names. For now, I at least put them next
#       to each other, so the duplication is more obvious.
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


def convert_table_comments_to_dict(tbl):

    comments_dict = {}
    if "comments" in tbl.meta:
        try:
            comments_str = "\n".join(tbl.meta["comments"])
            comments_dict = yaml.full_load(comments_str)
        except yaml.error.YAMLError:
            logger.warning("Couldn't convert <table>.meta['comments'] to dict")
            comments_dict = tbl.meta["comments"]
    elif "COMMENT" in tbl.meta:
        try:
            comments_dict = yaml.full_load("\n".join(tbl.meta["COMMENT"]))
        except yaml.error.YAMLError:
            logger.warning("Couldn't convert <table>.meta['COMMENT'] to dict")
            comments_dict = tbl.meta["COMMENT"]
    else:
        logger.debug("No comments in table")

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
            logger.warning("None of %s were found in %s", names, colnames)
    else:
        real_name = real_name[0]

    return real_name


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
    if is_bangkey(meta_dict[name]):
        raise ValueError(
            f"!-strings should be resolved upstream: {meta_dict[name]}")

    if isinstance(meta_dict[name], u.Quantity):
        unit = meta_dict[name].unit
    elif f"{name}_unit" in meta_dict:
        unit = meta_dict[f"{name}_unit"]
    else:
        unit = u.Unit(fallback_unit)

    return quantify(meta_dict[name], unit)


def quantify(item, unit, cmds=None):
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
        raise ValueError(f"Quantify cannot resolve {item}")
        # item = from_currsys(item, cmds)
    if isinstance(item, u.Quantity):
        return item.to(u.Unit(unit))

    if isinstance(item, (np.ndarray, list, tuple)) and np.size(item) > 1000:
        return item << u.Unit(unit)
    return item * u.Unit(unit)


def is_fits(filename) -> bool:
    # Using 'in ".fits"' to also catch ".fit", which exists sometimes...
    return (filename is not None and Path(filename).suffix.lower() in ".fits")


def get_fits_type(filename):
    with fits.open(filename) as hdulist:
        hdutype = "image"
        # pylint: disable=no-member
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
    logger.debug("%s_unit was not found in table.meta: %s. Default to: %s",
                 colname, tbl_name, default_unit)

    return u.Unit(default_unit)


def has_needed_keywords(header, suffix=""):
    """Check to see if the WCS keywords are in the header."""
    keys = {"CDELT1", "CRVAL1", "CRPIX1"}
    keys = {key + suffix for key in keys}
    keys.add("NAXIS1")
    return all(key in header.keys() for key in keys)


def stringify_dict(dic, ignore_types=(str, int, float, bool)):
    """Turn a dict entries into strings for addition to FITS headers."""
    for key, value in dic.items():
        if isinstance(value, ignore_types):
            yield key, value
        else:
            yield key, str(value)


def from_currsys(item, cmds=None):
    """Return the current value of a bang-string from ``rc.__currsys__``."""
    if isinstance(item, Table):
        tbl_dict = {col: item[col].data for col in item.colnames}
        tbl_dict = from_currsys(tbl_dict, cmds)
        item_meta = item.meta
        item = Table(data=list(tbl_dict.values()),
                     names=list(tbl_dict.keys()))
        item.meta = item_meta

    if isinstance(item, np.ndarray) and not isinstance(item, u.Quantity):
        item = np.array([from_currsys(x, cmds) for x in item])

    if isinstance(item, list):
        item = [from_currsys(x, cmds) for x in item]

    if isinstance(item, dict):
        for key in item:
            item[key] = from_currsys(item[key], cmds)

    if isinstance(item, str) and len(item) and item.startswith("!"):
        # if not isinstance(cmds, UserCommands)
        #     raise TypeError

        if not cmds:
            cmds = rc.__currsys__
            # raise ValueError(f"No cmds dict passed for resolving {item}")

        if item in cmds:
            item = cmds[item]
            if isinstance(item, str) and item.startswith("!"):
                item = from_currsys(item, cmds=cmds)
        else:
            raise ValueError(f"{item} was not found in rc.__currsys__")

    if isinstance(item, str):
        if item.lower() == "none":
            item = None
        try:
            item = float(item)
        except (TypeError, ValueError):
            pass

    return item


def from_rc_config(item):
    return from_currsys(item, rc.__config__)


def check_keys(input_dict: Union[Mapping, Iterable],
               required_keys: Set,
               action: str = "error",
               all_any: str = "all") -> bool:
    """
    Check to see if all/any of the required keys are present in a dict.

    .. versionchanged:: v0.8.0
        The `required_keys` parameter should now be a set.

    Parameters
    ----------
    input_dict : Union[Mapping, Iterable]
        The mapping to be checked.
    required_keys : Set
        Set containing the keys to look for.
    action : {"error", "warn", "warning"}, optional
        What to do in case the check does not pass. The default is "error".
    all_any : {"all", "any"}, optional
        Whether to check if "all" or "any" of the `required_keys` are present.
        The default is "all".

    Raises
    ------
    ValueError
        Raised when an invalid parameter was passed or when `action` was set to
        "error" (the default) and the `required_keys` were not found.

    Returns
    -------
    keys_present : bool
        ``True`` if check succeded, ``False`` otherwise.

    """
    # Checking for Set from collections.abc instead of builtin set to allow
    # for any duck typing (e.g. dict keys view or whatever)
    if not isinstance(required_keys, Set):
        logger.warning("required_keys should implement the Set protocol, "
                       "found %s instead.", type(required_keys))
        required_keys = set(required_keys)

    if all_any == "all":
        keys_present = required_keys.issubset(input_dict)
    elif all_any == "any":
        keys_present = not required_keys.isdisjoint(input_dict)
    else:
        raise ValueError("all_any must be either 'all' or 'any'")

    if not keys_present:
        missing = "', '".join(required_keys.difference(input_dict)) or "<none>"
        if "error" in action:
            raise ValueError(
                f"The keys '{missing}' are missing from input_dict.")
        if "warn" in action:
            logger.warning(
                "The keys '%s' are missing from input_dict.", missing)

    return keys_present


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
    # TODO: merge this functionality with the nested dict stuff in astar-utils
    text = ""
    for key, value in dic.items():
        if isinstance(value, dict):
            text += " " * indent + f"{str(key)}:\n"
            text += pretty_print_dict(value, indent=indent + 2)
        else:
            text += " " * indent + f"{str(key)}: {str(value)}\n"

    return text


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


def top_level_catch(func):
    """Catch any unhandled exceptions, log it including bug report."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            output = func(*args, **kwargs)
        except Exception as err:
            # FIXME: This try-except should not be necessary, but
            # logger.exception has an issue in some versions.
            try:
                bug_logger.exception(
                    "Unhandled exception occured, see log file for details.")
            except TypeError:
                bug_logger.error(
                    "Unhandled exception occured, see log file for details.")
                bug_logger.error("Couldn't log full exception stack.")
                bug_logger.error("Error message was: '%s'", err)
            log_bug_report(logging.ERROR)
            raise
        return output
    return wrapper


def update_logging(capture_warnings=True):
    """Reload logging configuration from ``rc.__logging_config__``."""
    # Need to access NestedMapping's internal dict here...
    dictConfig(rc.__logging_config__)
    logging.captureWarnings(capture_warnings)


def log_to_file(enable=True):
    """Enable or disable logging to file (convenience function)."""
    if enable:
        handlers = ["console", "file"]
    else:
        handlers = ["console"]

    rc.__logging_config__["loggers"]["astar"]["handlers"] = handlers
    update_logging()


def set_console_log_level(level="INFO"):
    """Set the level for the console handler (convenience function).

    This controls what is actually printed to the console by ScopeSim.
    Accepted values are: DEBUG, INFO (default), WARNING, ERROR and CRITICAL.
    """
    rc.__logging_config__["handlers"]["console"]["level"] = level
    update_logging()


def seq(start, stop, step=1):
    """Replacement for numpy.arange modelled after R's seq function.

    Returns an evenly spaced sequence from start to stop. stop is
    included if the difference between start and stop is an integer
    multiple of step.  From the documentation of numpy.range: "When
    using a non-integer step, such as 0.1, the results will often not
    be consistent." This replacement aims to avoid these
    inconsistencies.

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
