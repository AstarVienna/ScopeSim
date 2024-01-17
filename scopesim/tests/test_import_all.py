"""Test whether all Python files can be imported."""
import importlib.util
import time
from pathlib import Path

import scopesim


def test_import_all():
    """Test whether all Python files can be imported.

    This ensures that all files are included in the code coverage.
    """
    for name_path_scopesim in scopesim.__path__:
        path_scopesim = Path(name_path_scopesim)
        fns_python = path_scopesim.glob("**/*.py")
        for fn in fns_python:
            name_package = ".".join(fn.relative_to(path_scopesim.parent).with_suffix("").parts)
            time_1 = time.time()
            importlib.import_module(name_package)
            time_2 = time.time()
            time_delta = time_2 - time_1
            msg = f"{name_package} takes {time_delta} seconds to import"
            # Assert the package is quick to import, because this can highlight
            # broken packages. Requires running this test in isolation.
            assert time_delta < 0.2, msg
