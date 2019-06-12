import numpy as np
from astropy.io import fits
import scopesim as sim

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


# DETECTOR nested dictionary
# Ideally I would like to have "kwargs" : {"filename" : "MICADO_detectors.tbl"}
# so that MICADO_detectors.tbl is the one place where the detector description
# is stored. This however is more suited to storing the info in a database
DETECTOR_YAML = {"object": "detector",
                 "name": "test_detector",
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
                              "kwargs": {"value": 0.2}
                              }]
                 }


OBSERVATIONS_DICT = {"SIM_PIXEL_SCALE" : 0.004,     # because optical train still need this (stupidly)
                     "OBS_DIT" : 10,                # [sec]
                     "OBS_NDIT" : 1,                # Not yet implemented
                     "SIM_DETECTOR_YAML" : DETECTOR_YAML
                     }


def test_simplecado():

    src = sim.source.templates.empty_sky()

    cmd = sim.UserCommands(OBSERVATIONS_DICT)

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


def test_read_in_simplecado_package():
    sim.rc.__search_path__ = ["C:/Work/irdb/SimpleCADO/"]
    assert sim.utils.find_file("SimpleCADO.config")

    cmd = sim.UserCommands(filename="SimpleCADO.config")
    assert len(cmd.yaml_dicts) > 0

    opt = sim.OpticalTrain(cmds=cmd)
    assert opt.optics_manager.optical_elements[1].meta["object"] == "detector"

    src = sim.source.templates.empty_sky()
    opt.observe(src)
    hdu = opt.readout()

    assert type(hdu) == fits.HDUList

    if PLOTS:
        plt.imshow(hdu[5].data, norm=LogNorm())
        plt.show()


def test_download_simplecado_package():
    sim.rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"] = "C:/Work/simcado_downloads"

    dname = sim.rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"]
    sim.server.set_up_local_package_directory(dname, overwrite=True)

    sim.server.download_package("SimpleCADO")



