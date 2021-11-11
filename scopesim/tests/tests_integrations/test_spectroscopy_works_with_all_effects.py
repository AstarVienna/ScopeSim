import os
import pytest
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

from scopesim.optics import FOVManager, FieldOfView, ImagePlane
from scopesim.base_classes import PoorMansHeader
import scopesim.effects as efs
from scopesim.tests.mocks.py_objects import integr_spectroscopy_objects as iso

from scopesim import rc
MOCK_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        "../mocks/MICADO_SPEC/"))
rc.__search_path__.insert(0, MOCK_DIR)

PLOTS = False




################################################################################
# Everything needed to test the FOVManager in Spectroscopy mode


@pytest.fixture(scope="function")
def ap_list():
    return iso.mock_aperture_list()


@pytest.fixture(scope="function")
def ap_list_mixed():
    return iso.mock_aperture_list_mixed()\


@pytest.fixture(scope="function")
def ap_list_single():
    return iso.mock_aperture_list_single()


@pytest.fixture(scope="function")
def spt_list():
    return iso.mock_spectral_trace_list()


@pytest.fixture(scope="function")
def spt_list_shear():
    return iso.mock_spectral_trace_list_shear()


@pytest.fixture(scope="function")
def spt_list_single():
    return iso.mock_spectral_trace_list_single()


@pytest.fixture(scope="function")
def det_list():
    return iso.mock_detector_list()


@pytest.fixture(scope="function")
def config_yaml():
    return iso.mock_config_yaml()


@pytest.fixture(scope="function")
def point_source():
    return iso.mock_point_source_object()


@pytest.fixture(scope="function")
def ext_source():
    return iso.mock_extended_source_object()


@pytest.fixture(scope="function")
def gauss_psf():
    return iso.mock_gauss_psf()


@pytest.fixture(scope="function")
def shift_3d():
    return iso.mock_3d_shift()

################################################################################

@pytest.mark.skip(reason="Ignoring old Spectroscopy integration tests")
class TestSpectroscopyFOVs:
    @pytest.mark.usefixtures("ap_list", "spt_list", "det_list", "config_yaml",
                             "point_source", "ext_source", "gauss_psf")
    def test_basic_spectroscopy_mode(self, ap_list, spt_list, det_list,
                                     config_yaml, point_source, ext_source,
                                     gauss_psf):
        # Uses a basic SpectralTraceList with 4 traces and 2 apertures
        # 2 traces attached to each apertures
        config_yaml["wave_min"] = 1.0
        config_yaml["wave_max"] = 2.5

        src = point_source + ext_source

        fov_setup_effects = [ap_list, spt_list, det_list]
        fov_apply_effects = [gauss_psf]

        fov_mgr = FOVManager(effects=fov_setup_effects, **config_yaml)
        fovs = fov_mgr.fovs

        assert all([isinstance(fov, FieldOfView) for fov in fovs])

        implane = ImagePlane(det_list.image_plane_header)
        for fov in fovs:
            fov.extract_from(src)
            for effect in fov_apply_effects:
                fov = effect.apply_to(fov)

            if fov.hdu.data is None:
                fov.view()

            implane.add(fov.hdu, wcs_suffix="D")

        if PLOTS:
            plt.imshow(implane.data, origin="lower")
            plt.show()

    @pytest.mark.usefixtures("ap_list_mixed", "spt_list_shear", "det_list",
                             "config_yaml", "point_source", "ext_source",
                             "gauss_psf", "shift_3d")
    def test_spec_with_different_apertures(self, ap_list_mixed, spt_list_shear,
                                           det_list, config_yaml, point_source,
                                           ext_source, gauss_psf, shift_3d):
        # Uses a basic SpectralTraceList with 4 traces and 2 apertures
        # 2 traces attached to each apertures
        config_yaml["wave_min"] = 1.0
        config_yaml["wave_mid"] = 1.5
        config_yaml["wave_max"] = 2.5

        src = point_source + ext_source

        fov_setup_effects = [ap_list_mixed, spt_list_shear, det_list]
        fov_apply_effects = [gauss_psf, shift_3d]

        fov_mgr = FOVManager(effects=fov_setup_effects, **config_yaml)
        fovs = fov_mgr.fovs

        assert all([isinstance(fov, FieldOfView) for fov in fovs])

        implane = ImagePlane(det_list.image_plane_header)
        for fov in fovs:
            fov.extract_from(src)
            for effect in fov_apply_effects:
                fov = effect.apply_to(fov)

            if fov.hdu.data is None:
                fov.view()

            implane.add(fov.hdu, wcs_suffix="D")

        if PLOTS:
            plt.imshow(implane.data, origin="lower")
            plt.show()

    @pytest.mark.usefixtures("ap_list_single", "spt_list_single", "det_list",
                             "config_yaml", "point_source", "ext_source",
                             "gauss_psf", "shift_3d")
    def test_single_non_straight_traces(self, ap_list_single, spt_list_single,
                                        det_list, config_yaml, point_source,
                                        ext_source, gauss_psf):
        config_yaml["wave_min"] = 1.0
        config_yaml["wave_mid"] = 1.2
        config_yaml["wave_max"] = 2.5

        src = point_source + ext_source

        fov_setup_effects = [ap_list_single, spt_list_single, det_list]
        fov_apply_effects = [gauss_psf]

        fov_mgr = FOVManager(effects=fov_setup_effects, **config_yaml)
        fovs = fov_mgr.fovs

        assert all([isinstance(fov, FieldOfView) for fov in fovs])

        implane = ImagePlane(det_list.image_plane_header)
        for fov in fovs:
            fov.extract_from(src)
            for effect in fov_apply_effects:
                fov = effect.apply_to(fov)

            if fov.hdu.data is None:
                fov.view()

            implane.add(fov.hdu, wcs_suffix="D")

        if PLOTS:
            plt.imshow(implane.data, origin="lower")
            plt.show()

@pytest.mark.skip(reason="Ignoring old Spectroscopy integration tests")
class TestSpectroscopyMICADO:
    def test_initialises_spectral_trace_file(self):
        config = {"!SIM.spectral.wave_min": 1.45,
                  "!SIM.spectral.wave_mid": 1.9,
                  "!SIM.spectral.wave_max": 2.5,
                  "!INST.pixel_scale": 0.004,
                  "!INST.plate_scale": 0.26666667,
                  }
        for key in config:
            rc.__currsys__[key] = config[key]

        ap_mask = efs.ApertureMask(filename="SLIT_3000x50mas.dat")
        assert isinstance(ap_mask, efs.ApertureMask)

        spt = efs.SpectralTraceList(filename="TRACE_15arcsec.fits",
                                    wave_colname="lam", s_colname="xi",
                                    col_number_start=1)
        assert isinstance(spt, efs.SpectralTraceList)

        waves = spt.fov_grid(which="waveset")
        print("# waves:", len(waves))

        params = {"sky_header": ap_mask.header,
                  "det_header": None,
                  "pixel_scale": "!INST.pixel_scale",
                  "plate_scale": "!INST.plate_scale",
                  "wave_min": "!SIM.spectral.wave_min",
                  "wave_mid": "!SIM.spectral.wave_mid",
                  "wave_max": "!SIM.spectral.wave_max"}
        fovs = spt.fov_grid(which="edges", **params)
        print("# fovs:", len(fovs))
        assert isinstance(fovs[-1], (fits.Header, PoorMansHeader))
