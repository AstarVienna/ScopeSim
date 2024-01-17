"""Test whether all Python files can be imported."""
import importlib.util
from pathlib import Path

import scopesim


def test_import_all():
    """Test whether all Python files can be imported.

    This ensures that all files are included in the code coverage.
    """
    path_scopesim = Path(scopesim.__file__).parent
    paths_modules = path_scopesim.glob("**/*.py")
    for pm in paths_modules:
        spec = importlib.util.spec_from_file_location("module.name", pm)
        _ = importlib.util.module_from_spec(spec)
