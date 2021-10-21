import os
import yaml
from scopesim import rc

YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../yamls/"))
FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../files/"))
for NEW_PATH in [YAMLS_PATH, FILES_PATH]:
    if NEW_PATH not in rc.__search_path__:
        rc.__search_path__.insert(0, NEW_PATH)


def _atmo_yaml_dict():
    text = """
# ATMOSPHERE
object : atmosphere
name : armazones

properties :
    temperature : 0         # [-270..270] deg C
    pressure : 0.6          # [0..1] bar
    altitude : 3060
    longitude : -70.1918
    latitude : -24.5899
    airmass : 1.2
    humidity : 0.6
    pwv : 2.5
    pupil_angle : 30
    
effects :
-   name : super_psf
    class : GaussianDiffractionPSF
    kwargs :
        diameter : 39

-   name : atmo_dispersion
    description : atmospheric dispersion
    class : AtmosphericDispersion
    kwargs :
        wave_min : 1.9
        wave_min : 2.16
        wave_min : 2.4
        pixel_scale: 0.004
        
-   name : ignorable_effect
    class : Effect
    include : False
"""
    return yaml.full_load(text)


def _inst_yaml_dict():
    text = """
# INSTRUMENT OPTICS
object : instrument
name : micado_wide_field
inst_pkg_name : micado

properties :
    temperature : -190
    pixel_scale : 0.004
    

effects :
-   name : micado_surface_list
    class : SurfaceList
    kwargs :
        file_name : LIST_mirrors_MICADO_Wide.tbl

-   name : micado_adc
    class : AtmosphericDispersion
    kwargs :
        zenith_distance : 30
        reverse_shifts : True
        airmass : 1
        temperature : 0
        humidity : 0
        pressure : 0
        latitude : 0
        altitude : 0
        pupil_angle : 0
        
    """
    return yaml.full_load(text)


def _detector_yaml_dict():
    text = """
# Detector array
object : detector
name : micado_detector_array
inst_pkg_name : micado

properties :
    temperature : -190
    dark_current : 0.1
    image_plane_id : 0

effects :
-   name : detector_qe_curve
    class : TERCurve
    kwargs :
        filename : TER_blank.dat

-   name : micado_detector_geometry
    class : DetectorList
    kwargs:
        filename: LIST_detector_layout.dat
    """
    return yaml.full_load(text)


def _yaml_min_viable_scope():
    with open(os.path.join(YAMLS_PATH, "min_viable_sys.yaml")) as f:
        dicts = [dic for dic in yaml.full_load_all(f)]
    return dicts


def _usr_cmds_min_viable_scope():
    with open(os.path.join(FILES_PATH, "CMD_mvs_cmds.yaml")) as f:
        yaml_dicts = [dic for dic in yaml.full_load_all(f)]
    return yaml_dicts


def _yaml_unity_system():
    with open(os.path.join(YAMLS_PATH, "unity_sys.yaml")) as f:
        dicts = [dic for dic in yaml.full_load_all(f)]
    return dicts
