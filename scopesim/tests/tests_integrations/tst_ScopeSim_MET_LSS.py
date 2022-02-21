import pytest
import cProfile

import matplotlib.pyplot as plt
from scopesim.source import source_templates
from matplotlib.colors import LogNorm
import scopesim as sim
from scopesim import rc

rc.__currsys__['!SIM.file.local_packages_path'] = r"F:/Work/irdb"


def run_metis_lss():
    # src = sim.source.source_templates.empty_sky()
    spec = source_templates.ab_spectrum()
    src = sim.Source(x=[-1, 0, 1], y=[0, 0, 0],
                     ref=[0, 0, 0], weight=[1, 1, 1],
                     spectra=[spec])
    cmds = sim.UserCommands(use_instrument="METIS", set_modes=["lss_m"])
    metis = sim.OpticalTrain(cmds)
    metis["metis_psf_img"].include = False

    pr = cProfile.Profile()
    pr.enable()
    metis.observe(src)
    pr.disable()
    pr.print_stats(sort="cumulative")

    hdus = metis.readout()

    plt.subplot(122)
    plt.imshow(hdus[0][1].data, origin="lower", norm=LogNorm())
    plt.title("Detctor Plane (with noise)")
    plt.colorbar()

    plt.subplot(121)
    plt.imshow(metis.image_planes[0].data, origin="lower", norm=LogNorm())
    plt.title("Image Plane (noiseless)")
    plt.colorbar()
    plt.show()


run_metis_lss()


# class TestMetisLss:
#     run_metis_lss()
#
# def test_misc():
#     from astropy.io import fits
#     print(fits.info(r"F:\Work\irdb\METIS\TRACE_LSS_L.fits"))