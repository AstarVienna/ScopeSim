import os
import pytest
from astropy.table import Table
import telescopy.data.utils as tdb
import telescopy.default_keywords as dkeys


def test_data_database():
    assert "Dave" in tdb.data_database()

def test_get_local_pkgs_returns_Table_if_local_DB_file_exists():
    local_db_path = os.path.join(dkeys.PKG_DIR,
                                 dkeys.INST_PKG_LOCAL_PATH,
                                 dkeys.INST_PKG_LOCAL_DB_NAME)
    assert type(tdb.get_local_packages(local_db_path)) == Table

def test_get_local_pkgs_throws_error_if_local_DB_file_path_is_bogus():
    local_db_path = "bogus.txt"
    with pytest.raises(ValueError):
        tdb.get_local_packages(local_db_path)
