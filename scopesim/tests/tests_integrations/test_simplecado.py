import numpy as np
from astropy.io import fits
import scopesim as sim

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import scopesim.source.source_utils

PLOTS = False

# DETECTOR nested dictionary
# Ideally I would like to have "kwargs" : {"filename" : "MICADO_detectors.tbl"}
# so that MICADO_detectors.tbl is the one place where the detector description
# is stored. This however is more suited to storing the info in a database
DETECTOR_YAML = {"object": "detector",
                 "alias": "DET",
                 "name": "test_detector",
                 "properties": {"dit": "!OBS.dit"},
                 "effects": [{"name": "detector_array_list",
                              "description": "SimpleCADO detector array list",
                              "class": "DetectorList",
                              "kwargs": {"array_dict": {"id": [1],
                                                        "pixsize": [0.015],
                                                        "angle": [0.],
                                                        "gain": [1.0],
                                                        "x_cen": [0],
                                                        "y_cen":[0],
                                                        "xhw": [30.72],
                                                        "yhw": [30.72]},
                                         "x_cen_unit" : "mm",
                                         "y_cen_unit" : "mm",
                                         "xhw_unit" : "mm",
                                         "yhw_unit" : "mm",
                                         "pixsize_unit" : "mm",
                                         "angle_unit" : "deg",
                                         "gain_unit" : "electron/adu"
                                         }
                              },
                             {"name": "dark_current",
                              "description": "SimpleCADO dark current",
                              "class": "DarkCurrent",
                              # [e-/s] level of dark currentSimpleCADO.yaml
                              "kwargs": {"value": 0.2,
                                         "dit": "!OBS.dit",
                                         "ndit": "!OBS.ndit",}
                              }]
                 }

OBSERVATIONS_DICT = {"!OBS.ndit": 1,            # Not yet implemented
                     "!OBS.dit" : 10,           # [sec]
                     "!INST.pixel_scale": 0.004 # because optical train still need this (stupidly)
                    }


def test_simplecado():

    src = scopesim.source.source_utils.empty_sky()
    cmd = sim.commands.UserCommands(yamls=[DETECTOR_YAML],
                                    properties=OBSERVATIONS_DICT)

    opt = sim.OpticalTrain(cmd)
    opt.observe(src)
    hdu = opt.readout()

    # Finished - now just testing the output
    # No rogue photon torpedoes
    print(opt.image_plane.image)
    assert np.all(opt.image_plane.image) == 0

    # dark = 0.2 ct/s, DIT = 10s, NDIT = 1
    print(hdu[1].data)
    assert np.all(hdu[1].data == 2.0)


def read_in_simplecado_package():
    # ..todo: FIX THIS!!!!
    # only works on the local laptop
    scopesim.rc.__search_path__.insert(0, ["C:/Work/irdb/SimpleCADO/"])
    assert sim.utils.find_file("SimpleCADO.config")

    cmd = sim.UserCommands(filename="SimpleCADO.config")
    assert len(cmd.yaml_dicts) > 0

    opt = sim.OpticalTrain(cmds=cmd)
    assert opt.optics_manager.optical_elements[1].meta["object"] == "detector"

    src = scopesim.source.source_utils.empty_sky()
    opt.observe(src)
    hdu = opt.readout()

    assert type(hdu) == fits.HDUList

    if PLOTS:
        plt.imshow(hdu[5].data, norm=LogNorm())
        plt.show()
