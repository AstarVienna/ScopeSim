import pytest
import cProfile

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scopesim as sim
from scopesim import rc

rc.__currsys__['!SIM.file.local_packages_path'] = r"F:\Work\irdb"

#
# class TestMetisLss:
#     run_metis_lss()
#
# def test_misc():
#     from astropy.io import fits
#     print(fits.info(r"F:\Work\irdb\METIS\TRACE_LSS_L.fits"))

def run_metis_lss(pr):
    src = sim.source.source_templates.empty_sky()
    cmds = sim.UserCommands(use_instrument="METIS", set_modes=["lss_l"])
    metis = sim.OpticalTrain(cmds)
    metis["metis_psf_img"].include = False

    pr.enable()
    metis.observe(src)
    pr.disable()

    #hdus = metis.readout()

    # plt.subplot(122)
    # plt.imshow(hdus[0][1].data, origin="lower")
    # plt.subplot(121)
    # plt.imshow(metis.image_planes[0].data, origin="lower", norm=LogNorm())
    # plt.colorbar()
    # plt.show()

pr = cProfile.Profile()
run_metis_lss(pr)

# after your program ends
pr.print_stats(sort="cumulative")

