"""Test Source.shift by doing an integration test with LFAO."""
import pytest
import os
import shutil

import numpy
import scipy
from astropy import units as u
from numpy.testing import assert_approx_equal

import scopesim
from scopesim import rc
import scopesim_templates as sim_tp

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


if rc.__config__["!SIM.tests.run_integration_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring LFAO integration tests")

rc.__config__["!SIM.file.local_packages_path"] = "./lfoa_temp/"
rc.__config__["!SIM.file.use_cached_downloads"] = False

PKGS = {"LFOA": "telescopes/LFOA.zip"}

CLEAN_UP = False
PLOTS = False


def setup_module():
    """Download packages."""
    rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
    if not os.path.exists(rc_local_path):
        os.mkdir(rc_local_path)
        rc.__config__["!SIM.file.local_packages_path"] = os.path.abspath(
            rc_local_path)

    for pkg_name in PKGS:
        if not os.path.isdir(os.path.join(rc_local_path, pkg_name)) and \
                "irdb" not in rc_local_path:
            scopesim.download_package(PKGS[pkg_name])


def teardown_module():
    """Delete packages."""
    rc_local_path = rc.__config__["!SIM.file.local_packages_path"]
    if CLEAN_UP and "irdb" not in rc_local_path:
        shutil.rmtree(rc_local_path)


class TestShiftSource:
    def test_shift_lfao(self):
        # core_radius = 0.6 to ensure it fits the image after shifting
        src = sim_tp.basic.stars.cluster(mass=10000, distance=2000,
                                         core_radius=0.6,)

        lfoa = scopesim.OpticalTrain("LFOA")
        lfoa.observe(src)
        hdulists1 = lfoa.readout()

        data1 = hdulists1[0][1].data
        data = data1
        dmin, dmax, dmean, dmed, dstd = data.min(), data.max(), data.mean(), numpy.median(data), data.std()
        cm1y, cm1x = scipy.ndimage.center_of_mass(data1)
        print(cm1y, cm1x)

        if PLOTS:
            fig, ax = plt.subplots()
            im = ax.imshow(data, norm=LogNorm(vmin=dmed, vmax=dmed + 0.1 * dstd))
            fig.colorbar(im)
            fig.show()

        # Shift the cluster.
        dx = 10 * u.arcsec
        dy = 20 * u.arcsec
        src.shift(dx=dx, dy=dy)

        lfoa = scopesim.OpticalTrain("LFOA")
        lfoa.observe(src)
        hdulists2 = lfoa.readout()

        data2 = hdulists2[0][1].data
        data = data2
        dmin, dmax, dmean, dmed, dstd = data.min(), data.max(), data.mean(), numpy.median(data), data.std()
        if PLOTS:
            fig, ax = plt.subplots()
            im = ax.imshow(data, norm=LogNorm(vmin=dmed, vmax=dmed + 0.1 * dstd))
            # im = ax.imshow(data)
            fig.colorbar(im)
            fig.show()

        cm2y, cm2x = scipy.ndimage.center_of_mass(data2)
        print(cm2y, cm2x)

        # Compare the center of masses. Centers of mass. Centers of masses...
        dxm = cm2x - cm1x
        dym = cm2y - cm1y
        print(dxm, dym)

        # E.g. 10.0 == 10.1.
        assert_approx_equal(dx.value, dxm, 1)
        assert_approx_equal(dy.value, dxm, 1)
