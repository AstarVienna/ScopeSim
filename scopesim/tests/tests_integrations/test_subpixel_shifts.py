import pytest

from matplotlib import pyplot as plt
from astropy import units as u

from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim import OpticalTrain, UserCommands

PLOTS = False


@pytest.mark.usefixtures("protect_currsys")
def test_sub_pixels_integration():
    yaml = {"alias": "DET",
            "effects": [{"class": "DetectorWindow",
                         "kwargs": {"pixel_size": 0.015,
                                    "x": 0, "y": 0, "width": 0.165}
                         },
                        ]}
    properties = {"!SIM.sub_pixel.flag": True,
                  "!SIM.sub_pixel.fraction": 0.001,
                  "!INST.pixel_scale": 0.0015,
                  "!INST.plate_scale": 0.1}

    cmd = UserCommands(yamls=[yaml], properties=properties)
    opt = OpticalTrain(cmd)
    opt.cmds["!TEL.area"] = 1 * u.m**2

    src = so._vega_source(mag=15, x=0.001, y=0.)

    opt.observe(src)
    if PLOTS:
        plt.imshow(opt.image_planes[0].data, origin="lower")
        plt.show()
