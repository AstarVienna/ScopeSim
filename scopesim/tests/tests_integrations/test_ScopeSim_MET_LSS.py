import matplotlib.pyplot as plt
import scopesim as sim
from scopesim import rc

rc.__currsys__['!SIM.file.local_packages_path'] = r"F:\Work\irdb"


class TestMetisLss:
    def test_works(self):
        src = sim.source.source_templates.empty_sky()

        cmds = sim.UserCommands(use_instrument="METIS", set_modes=["lss_l"])
        metis = sim.OpticalTrain(cmds)

        metis["metis_psf_img"].include = False

        metis.observe(src)
        hdus = metis.readout()

        plt.subplot(121)
        plt.imshow(metis.image_planes[0].data, origin="lower")
        plt.subplot(122)
        plt.imshow(hdus[0][1].data, origin="lower")
        plt.show()


def test_misc():
    from astropy.io import fits
    print(fits.info(r"F:\Work\irdb\METIS\TRACE_LSS_L.fits"))