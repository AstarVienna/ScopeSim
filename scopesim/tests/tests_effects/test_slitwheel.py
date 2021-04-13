'''Tests for class SlitWheel'''
import os
import pytest

from scopesim import rc
from scopesim.effects import ApertureMask, SlitWheel

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
if FILES_PATH not in rc.__search_path__:
    rc.__search_path__ += [FILES_PATH]

@pytest.fixture(scope="class")
def swheel():
    '''Instantiate a SlitWheel'''
    return SlitWheel(slit_names=["A", "B"],
                     filename_format="MASK_slit_{}.dat",
                     current_slit="B")

class TestSlitWheel:
    def test_initialises_correctly(self, swheel):
        assert isinstance(swheel, SlitWheel)

    def test_has_aperture_masks(self, swheel):
        assert len(swheel.table) == 2

    def test_current_slit_is_aperture_mask(self, swheel):
        assert isinstance(swheel.current_slit, ApertureMask)
