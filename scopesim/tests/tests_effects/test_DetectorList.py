import os

from scopesim import rc
from scopesim.effects import DetectorList, ApertureMask

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
if MOCK_PATH not in rc.__search_path__:
    rc.__search_path__ += [MOCK_PATH]


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(DetectorList(), DetectorList)

    def test_initialises_with_filename(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0)
        assert isinstance(det_list, DetectorList)


class TestImagePlaneHeader:
    def test_header_is_sensical(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0)
        hdr_big = det_list.image_plane_header
        assert hdr_big["NAXIS1"] > 4096*3

    def test_header_fits_to_only_single_active_detector(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0,
                                active_detectors=[5])
        hdr_big = det_list.image_plane_header
        assert hdr_big["NAXIS1"] == 4096

    def test_header_fits_to_multiple_active_detectors(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0,
                                active_detectors=[1, 5])
        hdr_big = det_list.image_plane_header
        assert 4096 * 2 < hdr_big["NAXIS1"] < 4096 * 2 + 200
        assert 4096 * 2 < hdr_big["NAXIS2"] < 4096 * 2 + 200


class TestFovGrid:
    def test_returns_aperture_mask_object(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0)
        apm = det_list.fov_grid(pixel_scale=0.004)
        assert isinstance(apm, ApertureMask)

    def test_aperture_mask_header_covers_all_of_detector_header(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0)
        apm = det_list.fov_grid(pixel_scale=0.004)
        apm_hdr = apm.header
        det_hdr = det_list.image_plane_header
        assert apm_hdr["NAXIS1"] == det_hdr["NAXIS1"]
        assert apm_hdr["NAXIS2"] == det_hdr["NAXIS2"]


