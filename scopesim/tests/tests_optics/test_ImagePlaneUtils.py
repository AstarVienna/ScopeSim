import pytest
import numpy as np

from scopesim.optics import image_plane_utils as imp_utils


class TestSplitHeader:
    def test_is_the_header_split_into_the_right_amount_of_chunks(self):
        hdr = imp_utils.header_from_list_of_xy([-1.024, 1.024],
                                               [-1.024, 1.024], 0.002)
        hdrs = imp_utils.split_header(hdr, 128)
        area_sum = np.sum([hdr["NAXIS1"] * hdr["NAXIS2"] for hdr in hdrs])
        assert len(hdrs) == 64
        assert area_sum == hdr["NAXIS1"] * hdr["NAXIS2"]

    @pytest.mark.parametrize("x, y, pix", [(0.19, 0.2, 0.01),
                                           (2.19, 1.55, 0.01),
                                           (1.29, 1.2, 0.02),
                                           (2.55, 3.75, 1)])
    def test_final_header_is_smaller_for_odd_size(self, x, y, pix):
        hdr = imp_utils.header_from_list_of_xy([-1., x],
                                               [-1., y], pix)
        hdrs = imp_utils.split_header(hdr, 100)
        area_sum = np.sum([hdr["NAXIS1"] * hdr["NAXIS2"] for hdr in hdrs])
        assert area_sum == hdr["NAXIS1"] * hdr["NAXIS2"]

        # print([hdr["NAXIS1"] for hdr in hdrs], hdr["NAXIS1"])
        # print([hdr["NAXIS2"] for hdr in hdrs], hdr["NAXIS2"])
