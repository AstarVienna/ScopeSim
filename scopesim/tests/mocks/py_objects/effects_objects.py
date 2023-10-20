from unittest.mock import patch

from astropy import units as u
from astropy.table import Table

from scopesim import effects as efs
from scopesim.effects import AtmosphericDispersion
from scopesim.effects.effects_utils import make_effect
from scopesim.tests.mocks.py_objects.yaml_objects import _yaml_min_viable_scope

from . import FILES_PATH, MICADO_PATH


def _surf_list():
    kwargs = {"etendue": 5776 * u.m ** 2 * u.mas ** 2,
              "filename": str(MICADO_PATH / "LIST_mirrors_MICADO_Wide.tbl"),
              "name": "MICADO Mirror List"}
    with patch("scopesim.rc.__search_path__", [MICADO_PATH]):
        return efs.SurfaceList(**kwargs)


def _surf_list_empty():
    tbl = Table(names=["name", "outer", "inner", "angle",
                       "temperature", "action", "filename"],
                meta={"outer_unit": "m", "inner_unit": "m", "angle_unit": "deg",
                      "temperature_unit": "deg_C"})
    kwargs = {"etendue": 5776 * u.m ** 2 * u.mas ** 2,
              "table": tbl,
              "name": "Empty Surface List"}
    return efs.SurfaceList(**kwargs)


def _filter_surface(**kwargs):
    params = {"filename": str(MICADO_PATH / "TC_filter_Ks.dat"),
              "name": "filter",
              "action": "transmission",
              "outer": 0.1,
              "temperature": 0}
    params.update(kwargs)
    return efs.TERCurve(**params)


def _mvs_effects_list():
    effects_list = []
    for dic in _yaml_min_viable_scope():
        effects = dic["effects"]
        propeties = dic["properties"] if "properties" in dic else {}
        effects_list += [make_effect(eff, **propeties) for eff in effects]

    return effects_list


def _detector_list():
    kwargs = {"filename": str(FILES_PATH / "LIST_detector_layout.dat"),
              "image_plane_id": 0, "report": {}}
    return efs.DetectorList(**kwargs)


def _full_detector_list():
    kwargs = {"filename": str(FILES_PATH / "LIST_full_detector_layout.dat"),
              "image_plane_id": 0}
    return efs.DetectorList(**kwargs)


def _atmospheric_dispersion(**kwargs):
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


def _filter_tophat_curve():
    kwargs = {"name": "1-2um box-hat filter",
              "action": "transmission",
              "outer": 0.1,
              "temperature": 0,
              "wavelength_unit": "um",
              "array_dict": {"wavelength": [0.3, 0.99, 1., 2., 2.01, 3.0],
                             "transmission": [0, 0, 1, 1, 0, 0]}
              }
    return efs.FilterCurve(**kwargs)


def _const_psf():
    return efs.FieldConstantPSF(filename=str(FILES_PATH / "test_ConstPSF.fits"))


def _ncpa_psf():
    ncpa = efs.NonCommonPathAberration(pixel_scale=0.004)
    ncpa._total_wfe = 0.076
    return ncpa


def _img_aperture_mask(**kwargs):
    base_kwargs = {"array_dict": {"x": [-2, -1, 1, 2],
                                  "y": [-1, -2, 2, 1]},
                   "x_unit": "arcsec",
                   "y_unit": "arcsec"}
    base_kwargs.update(kwargs)
    apm = efs.ApertureMask(**base_kwargs)
    return apm
