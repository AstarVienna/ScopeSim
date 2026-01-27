import pytest
from pytest import approx

from matplotlib import pyplot as plt
import numpy as np

from scopesim.effects import RadialProfilePSF

PLOTS = False


class TestInit:
    def test_throws_error_for_no_profile(self):
        with pytest.raises(ValueError):
            RadialProfilePSF()

    def test_initalises_with_radial_profile(self):
        rppsf = RadialProfilePSF(r=[0, 3.5, 4.2],
                                 z=[1,   1,   0])
        assert isinstance(rppsf, RadialProfilePSF)

class TestApplyTo:
    def test_kernel_sums_to_unity(self):
        rppsf = RadialProfilePSF(r=[0, 3.5, 4.2],
                                 z=[1, 1, 0])
        kernel = rppsf.get_kernel(None)
        assert np.sum(kernel) == approx(1)

    def test_kernel_size_for_unit_pixel(self):
        rppsf = RadialProfilePSF(r=[0, 3.5, 4.2],
                                 z=[1, 1, 0])
        kernel = rppsf.get_kernel(None)
        assert kernel.shape[0] == 11 and kernel.shape[1] == 11

    def test_kernel_size_for_unit_mm(self):
        rppsf = RadialProfilePSF(r=[0, 3.5, 4.2],
                                 z=[1, 1, 0],
                                 unit="mm",
                                 pixel_scale=0.25,
                                 plate_scale=1)
        kernel = rppsf.get_kernel(None)
        assert kernel.shape[0] == int(np.ceil(4.2/0.25)*2+1)  # 35

        if PLOTS:
            plt.imshow(kernel)
            plt.show()

    def test_really_small_kernel_size(self):
        rppsf = RadialProfilePSF(r=[0, 0.9, 1.28, 1.5],
                                 z=[1, 1, 0.5, 0])
        kernel = rppsf.get_kernel(None)
        assert kernel.shape[0] == 5

        if PLOTS:
            plt.imshow(kernel)
            plt.show()
