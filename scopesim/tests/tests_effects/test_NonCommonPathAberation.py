import os

import pytest
from pytest import approx

from astropy import units as u

from scopesim.effects.psfs import NonCommonPathAberration
from scopesim.effects.psfs import strehl2gauss
from scopesim import rc

FILES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                         "../mocks/files"))
rc.__search_path__ = [FILES_DIR]


@pytest.fixture(scope="function")
def ncpa_kwargs():
    kwargs = {"pixel_scale": 0.004,
              "wfe_rms_unit": "nm",
              "array_dict": {"element": ["mirror", "entrance_window"],
                             "material": ["gold", "glass"],
                             "n_surfaces": [10, 18],
                             "wfe_rms": [20, 10]}}
    return kwargs


@pytest.mark.usefixtures("ncpa_kwargs")
class TestInit:
    def test_initialises_with_nothing_but_pixel_scale(self):
        ncpa = NonCommonPathAberration(pixel_scale=0.004)
        assert isinstance(ncpa, NonCommonPathAberration)

    def test_initialises_with_arrays(self, ncpa_kwargs):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        assert ncpa.table["wfe_rms"][0] == 20

    def test_initialises_with_file(self):
        kwargs = {"filename": "test_NCPAs_table.dat", "pixel_scale": 0.004}
        ncpa = NonCommonPathAberration(**kwargs)
        assert ncpa.table["wfe_rms"][0] == 20


@pytest.mark.usefixtures("ncpa_kwargs")
class TestGetKernel:
    def test_returns_total_wfe_in_units_of_table(self, ncpa_kwargs):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        ncpa.get_kernel([1.5, 2.5])
        assert ncpa.total_wfe.value == approx(67.082)

    def test_returns_kernel_with_proper_fwhm(self, ncpa_kwargs):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        ncpa.total_wfe = 1
        kernel = ncpa.get_kernel([1.5, 2.5]*u.um)

        import numpy as np
        print(np.max(kernel))

        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        plt.imshow(kernel, norm=LogNorm())
        plt.colorbar()
        # plt.show()


class TestStrehl2Gauss:
    def test_returns_delta_function_for_high_strehl(self):
        import numpy as np

        kernel = strehl2gauss(0.2222)
        print("")
        print(np.sum(kernel))
        print(np.max(kernel))

        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        plt.imshow(kernel)
        plt.colorbar()
        # plt.show()





