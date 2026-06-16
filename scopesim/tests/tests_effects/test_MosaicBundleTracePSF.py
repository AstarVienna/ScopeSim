import pytest
from pytest import approx

from matplotlib import pyplot as plt
import numpy as np

from scopesim.effects import MosaicBundleTracePSF, AnalyticalPSF

PLOTS = False


class TestInit:
    def test_initialised_with_nothing(self):
        assert isinstance(MosaicBundleTracePSF(), AnalyticalPSF)

    def test_initalises_with_radial_profile(self):
        mbt_psf = MosaicBundleTracePSF(r=[0, 3.5, 4.2],
                                     z=[1,   1,   0])
        assert isinstance(mbt_psf, MosaicBundleTracePSF)


class TestGetKernel:
    def test_returned_kernel_is_zipped_version_of_single_kernel(self):
        mbt_psf = MosaicBundleTracePSF(n_traces_per_bundle=2,
                                       trace_spacing=2.1,
                                       trace_flux_weights=[0.5, 2, 1])
        kernel = mbt_psf.get_kernel(None)

        if PLOTS:
            plt.imshow(kernel)
            plt.show()

        ws = mbt_psf.meta["trace_flux_weights"]
        assert np.sum(kernel) == np.sum(ws) * np.sum(mbt_psf.single_kernel)
        assert kernel.shape[0] > ws * kernel.shape[0]
