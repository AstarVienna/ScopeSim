from scopesim.effects import ApertureMask


def _basic_aperture():
    kwargs = {"array_dict": {"x": [-2, 2, 2, -2],
                             "y": [-0.1, -0.1, 0.1, 0.1]},
              "x_unit": "arcsec",
              "y_unit": "arcsec",
              "pixel_scale": 0.1}
    apm = ApertureMask(**kwargs)

    return apm


def _basic_micado_slit_aperture():
    kwargs = {"array_dict": {"x": [-1.5, 1.5, 1.5, -1.5],
                             "y": [-0.025, -0.025, 0.025, 0.025]},
              "x_unit": "arcsec",
              "y_unit": "arcsec",
              "pixel_scale": 0.004,
              "no_mask": False}
    apm = ApertureMask(**kwargs)

    return apm
