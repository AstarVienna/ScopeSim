from scopesim.detector import Detector
from scopesim.tests.mocks.py_objects.header_objects import _basic_dtcr_header


def _basic_detector():
    dtcr = Detector(_basic_dtcr_header())
    return dtcr
