import os
import pytest
from pytest import approx

from astropy import units as u

import scopesim as sim
from scopesim.optics.effects import DetectorList

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__ += [MOCK_PATH]


class TestDetectorInit:
    def test_initialises_with_nothing(self):
        assert isinstance(DetectorList(), DetectorList)

    def test_initialises_with_filename(self):
        det_list = DetectorList(filename="FPA_array_layout.dat")
        assert isinstance(det_list, DetectorList)


class TestDetectorImagePlaneHeader:
    def test_header_is_sensical(self):
        det_list = DetectorList(filename="FPA_array_layout.dat")
        hdr_big = det_list.image_plane_header
        assert hdr_big["NAXIS1"] > 4096*3
