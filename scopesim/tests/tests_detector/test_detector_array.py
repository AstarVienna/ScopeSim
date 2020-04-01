import pytest

import numpy as np
from astropy.io import fits

from scopesim.optics.image_plane import ImagePlane
from scopesim.detector import DetectorArray
from scopesim.tests.mocks.py_objects import effects_objects as efs_objs


@pytest.fixture(scope="function")
def detector_list_effect():
    return efs_objs._detector_list()


@pytest.fixture(scope="function")
def image_plane():
    dtcr_list = efs_objs._detector_list()
    implane = ImagePlane(header=dtcr_list.image_plane_header)
    implane.hdu.data = np.zeros((implane.header["NAXIS2"],
                                 implane.header["NAXIS1"]))
    return implane


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(DetectorArray(), DetectorArray)

    def test_initialisation_parameters_stored_in_meta(self):
        dtcr_arr = DetectorArray(random_parameter=42)
        assert dtcr_arr.meta["random_parameter"] == 42


@pytest.mark.usefixtures("image_plane", "detector_list_effect")
class TestReadout:
    def test_returns_hdu_for_empty_effects_list(self, image_plane,
                                                detector_list_effect):
        dtcr_arr = DetectorArray(detector_list_effect)
        hdu = dtcr_arr.readout([image_plane])

        assert isinstance(hdu, fits.HDUList)

    def test_hdu_data_is_lots_of_zeros_for_empty_input(self, image_plane,
                                                       detector_list_effect):
        dtcr_arr = DetectorArray(detector_list_effect)
        hdu = dtcr_arr.readout([image_plane])

        assert np.all(hdu[1].data == 0)
        assert hdu[1].shape[0] == detector_list_effect.table["x_len"]
