from scopesim.detector import Detector
from scopesim.tests.mocks.py_objects.header_objects import _basic_dtcr_header


def _basic_detector(width=16, pix_size=0.01):
    dtcr = Detector(_basic_dtcr_header(width, pix_size))
    return dtcr
