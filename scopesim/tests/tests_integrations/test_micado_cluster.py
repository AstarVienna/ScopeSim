"""Simulate a large cluster with ScopeSIM.

Test for https://github.com/AstarVienna/ScopeSim/issues/111

This code works fine in ScopeSIM 0.1.4, but breaks in ScopeSIM 0.4.0.

Error:

Traceback (most recent call last):
  File "scopesim_large_cluster_test.py", line 103, in <module>
    opt.observe(src)
  File "../site-packages/scopesim/optics/optical_train.py", line 188, in observe
    fov.view(hdu_type)
  File "../site-packages/scopesim/optics/fov.py", line 144, in view
    self.hdu = self.make_image_hdu(use_photlam=use_photlam)
  File "../site-packages/scopesim/optics/fov.py", line 338, in make_image_hdu
    canvas_image_hdu.data[y, x] += f * weight
IndexError: index 4096 is out of bounds for axis 0 with size 4096

Note that the error goes away the detector is reduced in size significantly.
E.g. by choosing xhw and yhw of about 0.72.
"""

import copy

import scopesim
import scopesim_templates
from collections import OrderedDict

import astropy.units as u

def test_cluster():
    scopesim_dictionaries = [
        {
            "object": "observation",
            "alias": "OBS",
            "name": "MICADO_configuration",
            "description": "parameters needed for a MICADO simulation",
            "properties": {"dit": 1.0, "ndit": 1},
        },
        {
            "class": "DetectorArray",
            "description": "A detector array.",
            "alias": "DET",
            "effects": [
                {
                    "class": "DetectorList",
                    "description": "A single Detector.",
                    "kwargs": {
                        "array_dict": OrderedDict(
                            [
                                ("id", [""]),
                                ("angle", [0.0]),
                                ("gain", [1.0]),
                                ("pixsize", [0.015]),
                                ("x_cen", [0.0]),
                                ("y_cen", [0.0]),
                                # Detector size must be large!
                                ("xhw", [30.72]),
                                ("yhw", [30.72]),
                            ]
                        ),
                        "angle_unit": "degree",
                        "gain_unit": "e/adu",
                        "pixsize_unit": "mm",
                        "x_cen_unit": "mm",
                        "y_cen_unit": "mm",
                        "xhw_unit": "mm",
                        "yhw_unit": "mm",
                    },
                    "name": "detector_array_list",
                },
            ],
            "properties": {
                "name": "testdetector",
                "image_plane_id": 0,
                "temperature": -230.0,
                "temperature_unit": "C",
                "dit": "!OBS.dit",
                "ndit": "!OBS.ndit",
            },
        },
    ]

    observations_dict = {
        "!INST.pixel_scale": 0.004,
        "!INST.plate_scale": 0.2666667,
    }

    src = scopesim_templates.micado.cluster()
    dx, dy = 0 * u.arcsec, 0 * u.arcsec
    # Source must be shifted!
    src.shift(dx, dy)
    cmd = scopesim.UserCommands(
        yamls=copy.copy(scopesim_dictionaries),
        properties=copy.copy(observations_dict),
    )
    opt = scopesim.OpticalTrain(cmd)
    opt.observe(src)
    # This should not raise an exception.
