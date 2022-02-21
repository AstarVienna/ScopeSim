import os

import pytest
from pytest import approx

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u

from scopesim.effects.psfs import NonCommonPathAberration
from scopesim.effects.psf_utils import strehl2sigma, sigma2gauss, wfe2gauss, wfe2strehl
from scopesim.optics import FieldOfView, ImagePlane
from scopesim import rc

from scopesim.tests.mocks.py_objects.source_objects import _single_table_source
from scopesim.tests.mocks.py_objects.header_objects import \
    _fov_header, _implane_header


FILES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                         "../mocks/files"))
if rc.__search_path__[0] != FILES_DIR:
    rc.__search_path__.insert(0, FILES_DIR)

PLOTS = False


@pytest.fixture(scope="function")
def fov_Ks():
    _src = _single_table_source()
    _fov = FieldOfView(header=_fov_header(), waverange=(1.9, 2.4), area=1*u.m**2)
    _fov.extract_from(_src)
    _fov.view()
    return _fov


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


@pytest.mark.usefixtures("ncpa_kwargs", "fov_Ks")
class TestGetKernel:
    def test_returns_total_wfe_in_units_of_table(self, ncpa_kwargs, fov_Ks):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        assert ncpa.total_wfe.value == approx(76.158, rel=1e-5)

    @pytest.mark.parametrize("waves, strehl", [((1.1, 1.3),  0.855),
                                               ((1.45, 1.8), 0.923),
                                               ((1.9, 2.4),  0.953)])
    def test_returns_kernel_with_expected_central_peak_value(self, fov_Ks,
                                                             ncpa_kwargs,
                                                             waves, strehl):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        fov_Ks.meta["wave_min"] = waves[0]
        fov_Ks.meta["wave_max"] = waves[1]
        kernel = ncpa.get_kernel(fov_Ks)

        assert np.abs(np.max(kernel) / strehl - 1) < 0.01
        assert kernel.shape == (3, 3)
        assert np.sum(kernel) == approx(1)

        if PLOTS:
            plt.imshow(kernel, norm=LogNorm(), vmax=1, vmin=1e-4)
            plt.colorbar()
            plt.show()


@pytest.mark.usefixtures("ncpa_kwargs", "fov_Ks")
class TestApplyTo:
    def test_convolves_kernel_with_fov_image(self, ncpa_kwargs, fov_Ks):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        pre_max_flux = np.max(fov_Ks.data)
        fov_Ks = ncpa.apply_to(fov_Ks)
        post_max_flux = np.max(fov_Ks.data)

        assert post_max_flux/pre_max_flux == approx(0.954, rel=0.002)

        if PLOTS:
            plt.imshow(fov_Ks.image[40:60, 40:60], norm=LogNorm())
            plt.show()

    def test_ignores_classes_other_than_fov(self, ncpa_kwargs):
        implane = ImagePlane(_implane_header())
        implane.hdu.data[1, 1] = 1
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        new_implane = ncpa.apply_to(implane)

        assert np.all(new_implane.data == implane.data)

        if PLOTS:
            plt.imshow(new_implane.image[:5,:5])
            plt.show()


@pytest.mark.usefixtures("ncpa_kwargs", "fov_Ks")
class TestFovGrid:
    def test_returns_currsys_edge_waves_for_no_input(self, ncpa_kwargs, fov_Ks):
        ncpa = NonCommonPathAberration(**ncpa_kwargs)
        waves = ncpa.fov_grid()
        wave_min = rc.__currsys__["!SIM.spectral.wave_min"]
        wave_max = rc.__currsys__["!SIM.spectral.wave_max"]
        assert waves[0].to(u.um).value == approx(wave_min)
        assert waves[-1].to(u.um).value == approx(wave_max)


################################################################################

class TestFunctionStrehl2Gauss:
    def test_relationship_between_sigma_strehl_amplitude(self):
        # test that the central pixel is equal to the strehl ratio needed

        wave = np.arange(0.3, 4, 0.1)
        srs = wfe2strehl(0.076, wave)
        sigs = strehl2sigma(srs)
        kernels = np.array([sigma2gauss(sig) for sig in sigs])
        amplis = np.array([np.max(kernel) for kernel in kernels])
        ampli_srs = amplis / srs

        assert np.all(0.96 < ampli_srs) and np.all(ampli_srs < 1.04)

        if PLOTS:
            plt.plot(amplis, srs, c="b")
            plt.show()

    def test_returns_expected_micado_strehl_values(self):
        """
        95.4% strehl at 2.2um and 91.9% at 1.65um and 86.4% at 1.25um
        """
        strehls = np.array([0.864, 0.919, 0.954])
        kernels = [wfe2gauss(0.076, wave) for wave in [1.25, 1.65, 2.2]]
        max_kernels = np.array([np.max(kernel) for kernel in kernels])
        ratios = max_kernels / strehls

        assert np.all(np.abs(ratios-1) < 0.01)
        assert np.all([k.shape == (3, 3) for k in kernels])
        assert np.all([np.sum(k) == approx(1) for k in kernels])
