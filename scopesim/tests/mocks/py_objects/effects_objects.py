import os

from astropy import units as u
from astropy.table import Table

from scopesim import effects as efs
from scopesim.effects.effects_utils import make_effect
from scopesim.tests.mocks.py_objects.yaml_objects import _yaml_min_viable_scope

from scopesim import rc

YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../files/"))
if YAMLS_PATH not in rc.__search_path__:
    rc.__search_path__ += [YAMLS_PATH]


def _surf_list():
    kwargs = {"etendue": 5776 * u.m ** 2 * u.mas ** 2,
              "filename": "LIST_mirrors_MICADO_Wide.tbl",
              "name": "MICADO Mirror List"}
    return efs.SurfaceList(**kwargs)


def _surf_list_empty():
    tbl = Table(names=["Name", "Outer", "Inner", "Angle",
                       "Temp", "Action", "Filename"],
                meta={"outer_unit": "m", "inner_unit": "m", "angle_unit": "deg",
                      "temp_unit": "deg_C"})
    kwargs = {"etendue": 5776 * u.m ** 2 * u.mas ** 2,
              "table": tbl,
              "name": "Empty Surface List"}
    return efs.SurfaceList(**kwargs)


def _filter_surface():
    kwargs = {"filename": "TC_filter_Ks.dat",
              "name": "filter",
              "action": "transmission",
              "outer": 0.1,
              "temp": 0}
    return efs.TERCurve(**kwargs)


def _mvs_effects_list():
    effects_list = []
    for dic in _yaml_min_viable_scope():
        effects = dic["effects"]
        effects_list += [make_effect(eff) for eff in effects]

    return effects_list


def _detector_list():
    kwargs = {"filename": "LIST_detector_layout.dat"}
    return efs.DetectorList(**kwargs)


def _atmospheric_dispersion(**kwargs):
    from scopesim.effects import AtmosphericDispersion
    atmo_params = {"airmass": 1.14,     # in deg
                   "temperature": 7,    # in degC
                   "humidity": 1,       # in %
                   "pressure": 0.755,   # in mbar
                   "latitude": -26.5,   # in deg
                   "altitude": 2400,    # in m
                   "wave_min": 0.5,     # in um
                   "wave_mid": 1.5,
                   "wave_max": 2.5,
                   "pixel_scale": 0.004, # in arcsec
                   "pupil_angle": 0,     # in deg
                   "sub_pixel_fraction": 1,
                   "num_steps": 1000}
    atmo_params.update(kwargs)
    return AtmosphericDispersion(**atmo_params)
