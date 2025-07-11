import pytest
from unittest.mock import patch

from astropy.table import Table
from astropy import units as u

from scopesim.effects import DetectorList, DetectorWindow, DetectorList3D


@pytest.fixture(scope="class")
def patch_mock_path_micado(mock_path_micado):
    with patch("scopesim.rc.__search_path__", [mock_path_micado]):
        yield


class TestDetectorListInit:
    def test_initialises_with_nothing(self):
        assert isinstance(DetectorList(), DetectorList)

    @pytest.mark.usefixtures("patch_mock_path_micado")
    def test_initialises_with_filename(self):
        det_list = DetectorList(filename="FPA_array_layout.dat",
                                image_plane_id=0)
        _ = det_list.detector_headers()[0]
        assert isinstance(det_list, DetectorList)
        assert "x_size" in det_list.table.colnames
        assert det_list.table["x_size"][0] == 61.44

    def test_initialised_with_table(self):
        x, y, width, height, angle, gain, pixel_size = 0, 0, 10, 10, 0, 1, 0.1
        tbl = Table(data=[[0], [x], [y], [width], [height],
                          [angle], [gain], [pixel_size]],
                    names=["id", "x_cen", "y_cen", "x_size", "y_size",
                           "angle", "gain", "pixel_size"])
        params = {"x_cen_unit": "mm", "y_cen_unit": "mm",
                  "x_size_unit": "mm", "y_size_unit": "mm"}
        det_list = DetectorList(table=tbl, **params)
        hdr = det_list.detector_headers()[0]
        assert hdr["NAXIS1"] == 100

    def test_with_old_column_names(self):
        x, y, width, height, angle, gain, pixel_size = 0, 0, 5, 5, 0, 1, 0.1
        tbl = Table(data=[[0], [x], [y], [width], [height],
                          [angle], [gain], [pixel_size]],
                    names=["id", "x_cen", "y_cen", "xhw", "yhw",
                           "angle", "gain", "pixsize"])
        params = {"x_cen_unit": "mm", "y_cen_unit": "mm",
                  "x_size_unit": "mm", "y_size_unit": "mm"}
        det_list = DetectorList(table=tbl, **params)
        hdr = det_list.detector_headers()[0]
        assert hdr["NAXIS1"] == 100


@pytest.mark.usefixtures("patch_mock_path_micado")
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

    def test_happy_using_x_size_unit_in_pixels(self):
        det_list = DetectorList(array_dict={
            "id": [0, 1],
            "x_cen": [-0.2, 0.2],     # mm
            "y_cen": [0, 0],
            "x_size": [32, 32],       # pixel
            "y_size": [32, 32],
            "pixsize": [0.01, 0.01],  # mm/pixel
            "angle": [0, 0],
            "gain": [1, 1]},
            x_size_unit="pixel",
            y_size_unit="pixel",
            image_plane_id=0)

        hdr_big = det_list.image_plane_header
        assert hdr_big["NAXIS1"] == (16 + 20) * 2
        assert hdr_big["NAXIS2"] == 32

    def test_happy_using_x_size_unit_in_mm(self):
        det_list = DetectorList(array_dict={
            "id": [0, 1],
            "x_cen": [-0.2, 0.2],     # mm
            "y_cen": [0, 0],
            "x_size": [0.32, 0.32],   # mm
            "y_size": [0.32, 0.32],
            "pixsize": [0.01, 0.01],  # mm/pixel
            "angle": [0, 0],
            "gain": [1, 1]},
            x_size_unit="mm",
            y_size_unit="mm",
            image_plane_id=0)

        hdr_big = det_list.image_plane_header
        assert hdr_big["NAXIS1"] == (16 + 20) * 2
        assert hdr_big["NAXIS2"] == 32


class TestDetecotrHeaders:
    def test_happy_using_x_size_unit_in_pixels(self):
        det_list = DetectorList(array_dict={
            "id": [0, 1],
            "x_cen": [-0.2, 0.2],     # mm
            "y_cen": [0, 0],
            "x_size": [32, 32],       # pixel
            "y_size": [32, 32],
            "pixsize": [0.01, 0.01],  # mm/pixel
            "angle": [0, 0],
            "gain": [1, 1]},
            x_size_unit="pixel",
            y_size_unit="pixel",
            x_cen_unit="mm",
            y_cen_unit="mm",
            image_plane_id=0)

        det_hdrs = det_list.detector_headers()
        for hdr in det_hdrs:
            assert hdr["NAXIS1"] == 32
            assert hdr["NAXIS2"] == 32

    def test_happy_using_x_size_unit_in_mm(self):
        det_list = DetectorList(array_dict={
            "id": [0, 1],
            "x_cen": [-0.2, 0.2],     # mm
            "y_cen": [0, 0],
            "x_size": [0.32, 0.32],   # mm
            "y_size": [0.32, 0.32],
            "pixsize": [0.01, 0.01],  # mm/pixel
            "angle": [0, 0],
            "gain": [1, 1]},
            x_size_unit="mm",
            y_size_unit="mm",
            x_cen_unit="mm",
            y_cen_unit="mm",
            image_plane_id=0)

        det_hdrs = det_list.detector_headers()
        for hdr in det_hdrs:
            assert hdr["NAXIS1"] == 32
            assert hdr["NAXIS2"] == 32


class TestDetectorWindowInit:
    def test_throws_error_when_initialises_with_nothing(self):
        with pytest.raises(TypeError):
            DetectorWindow()

    def test_initialises_with_correct_parameters(self):
        det_window = DetectorWindow(pixel_size=0.1, x=0, y=0, width=10)
        assert det_window.data["x_cen"][0] == 0
        assert det_window.data["y_size"][0] == 10
        assert det_window.data["pixel_size"][0] == 0.1

    def test_recognised_as_detector_list(self):
        det_window = DetectorWindow(pixel_size=0.1, x=0, y=0, width=10)
        assert isinstance(det_window, DetectorWindow)
        assert isinstance(det_window, DetectorList)

    def test_initialises_with_correct_new_col_name_x_size(self):
        det_window = DetectorWindow(pixel_size=0.1, x=0, y=0, width=10)
        assert "xhw" not in det_window.data.colnames
        assert len(det_window.detector_headers()) == 1
        assert det_window.data["x_size"][0] == 10

    def test_can_use_bang_strings_to_define_size(self):
        patched = {"!DET.width": 4.2,
                   "!DET.pixel_size": 0.1}
        with patch.dict("scopesim.rc.__currsys__", patched):
            det = DetectorWindow(pixel_size="!DET.pixel_size", x=0, y=0,
                                 width="!DET.width")
            assert det.detector_headers()[0]["NAXIS1"] == 42.

        patched["!DET.width"] = 900.1
        with patch.dict("scopesim.rc.__currsys__", patched):
            assert det.image_plane_header["NAXIS1"] == 9001

    def test_can_define_everything_in_pixels(self):
        patched = {"!DET.width": 42}
        with patch.dict("scopesim.rc.__currsys__", patched):
            det = DetectorWindow(pixel_size=0.1, x=0, y=0,
                                 width="!DET.width", units="pixel")
            assert det.image_plane_header["NAXIS1"] == 42
            assert det.detector_headers()[0]["NAXIS1"] == 42


@pytest.fixture(name="det_list")
def basic_detectorlist3d():
    det_list = DetectorList3D(array_dict={
        "id": [0],
        "x_cen": [0],
        "y_cen": [0],
        "z_cen": [0],
        "x_size": [15],
        "y_size": [25],
        "z_size": [5],
        "pixel_size": [0.01],
        "angle": [0],
        "gain": [1]},
        x_size_unit="pixel",
        y_size_unit="pixel",
        z_size_unit="pixel",
        image_plane_id=0,
        pixel_scale=0.5,
    )
    return det_list


@pytest.mark.usefixtures("patch_all_mock_paths")
class TestDetectorList3D:
    def test_image_plane_header_is_3d(self, det_list):
        impl_hdr = det_list.image_plane_header

        assert impl_hdr["NAXIS"] == 3
        assert impl_hdr["NAXIS1"] == 15
        assert impl_hdr["NAXIS2"] == 25
        assert impl_hdr["NAXIS3"] == 5

    def test_pixel_scale_arcsec(self, det_list):
        equiv = det_list.pixel_scale_arcsec
        assert (1 * u.pixel).to(u.arcsec, equiv) == .5 * u.arcsec

    def test_pixel_scale_mm(self, det_list):
        equiv = det_list.pixel_scale_mm
        assert (1 * u.pixel).to(u.mm, equiv) == .01 * u.mm
