import cProfile
import scopesim


def profile_hawki_integration():
    PKGS = {"Paranal": "locations/Paranal.zip",
            "VLT": "telescopes/VLT.zip",
            "HAWKI": "instruments/HAWKI.zip"}
    for pkg_name in PKGS:
        scopesim.download_package(PKGS[pkg_name])

    cmd = scopesim.UserCommands(use_instrument="HAWKI",
                                properties={"!OBS.dit": 360,
                                            "!OBS.ndit": 10})
    cmd.ignore_effects += ["detector_linearity"]

    opt = scopesim.OpticalTrain(cmd)
    src = scopesim.source.source_templates.star_field(10000, 5, 15, 440)

    # ETC gives 2700 e-/DIT for a 1s DET at airmass=1.2, pwv=2.5
    # background should therefore be ~ 8.300.000
    opt.observe(src)
    hdu = opt.readout()[0]


cProfile.run('profile_hawki_integration()')