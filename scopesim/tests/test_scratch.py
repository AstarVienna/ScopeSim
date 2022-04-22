import scopesim as sim
from astropy import units as u

import os

def test_lms():
    PKGS_DIR = os.path.abspath("../../../irdb/")
    sim.rc.__config__["!SIM.file.local_packages_path"] = PKGS_DIR

    cmd = sim.UserCommands(use_instrument="METIS", set_modes=['lms'])
    metis = sim.OpticalTrain(cmd)

    src = sim.source.source_templates.star(flux=0.25*u.Jy)

    metis.observe(src, update=True)