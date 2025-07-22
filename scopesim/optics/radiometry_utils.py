from copy import deepcopy

from .surface import SpectralSurface
from ..utils import real_colname, get_logger


logger = get_logger(__name__)


def combine_emissions(tbl, surfaces, row_indexes, etendue, use_area=False):
    """
    Combine thermal emission from a series of surfaces.

    The function traces thermal emission through an optical system, taking
    into account the finite reflectivities/transmissivities and emissivities
    of the surfaces. The function assumes that etendue is conserved through
    the system, i.e. surfaces are neither over- nor undersized.

    Parameters
    ----------
    tbl : astropy Table
        Required columns are `name` and `action` (reflection or transmission)
    surfaces: OrderedDict of SpectralSurface
        Keys are the names from tbl, values are of type `SpectralSurface`
    row_indexes : list of int
        Rows of tbl (i.e. surfaces) to combine
    etendue, use_area : not needed  (TODO: remove)

    Returns
    -------
    SourceSpectrum

    """
    if len(tbl) == 0:
        return None

    r_name = real_colname("name", tbl.colnames)
    r_action = real_colname("action", tbl.colnames)

    emission = None
    for row_num in row_indexes:
        row = tbl[row_num]
        surf = surfaces[row[r_name]]
        action_attr = row[r_action]

        if isinstance(surf, SpectralSurface):
            surf_throughput = getattr(surf, action_attr)
            if emission is not None:
                emission = emission * surf_throughput

            if emission is None:
                emission = deepcopy(surf.emission)
            else:
                emission = emission + surf.emission

    return emission


def combine_throughputs(tbl, surfaces, rows_indexes):
    if len(tbl) == 0:
        return None

    r_name = real_colname("name", tbl.colnames)
    r_action = real_colname("action", tbl.colnames)

    throughput = None
    for ii, row_num in enumerate(rows_indexes):

        row = tbl[row_num]
        surf = surfaces[row[r_name]]
        action_attr = row[r_action]
        if action_attr == "":
            raise ValueError(f"No action in surf.meta: {surf.meta}")

        if isinstance(surf, SpectralSurface):
            surf_throughput = getattr(surf, action_attr)

            if ii == 0:
                throughput = deepcopy(surf_throughput)
            else:
                throughput = throughput * surf_throughput

    return throughput
