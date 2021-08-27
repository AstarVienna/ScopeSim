import pytest
from astropy import units as u

from scopesim.tests.mocks.py_objects import header_objects as ho
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim.optics.fov2 import FieldOfView


class TestExtractFrom:
    def test_extract_point_sources_from_table(self):
        src = so._table_source()
        src.fields[0]["x"] = [-15,-5,0,0]
        src.fields[0]["y"] = [0,0,5,15]

        hdr = ho._fov_header()              # 20x20" @ 0.2" --> [-10, 10]"
        wav = [1.9, 2.1] * u.um

        fov = FieldOfView(hdr, wav)
        fov.extract_from(src)

        assert len(fov.fields[0]) == 2
