# -*- coding: utf-8 -*-
"""Global fixtures for pytest."""

import logging
from pathlib import Path

import pytest
from unittest.mock import patch

import scopesim as sim
from astar_utils import UniqueList


MOCK_DIR = Path(__file__).parent / "mocks"

sim.rc.__currsys__["!SIM.file.error_on_missing_file"] = True


@pytest.fixture(scope="function", autouse=True)
def configure_logging():
    """This disables the handlers so the logs reach pytests' caplog.

    This fixture should be on the function level, because some functions
    might call `update_logging()` and therefore undo this fixture for
    subsequent tests. E.g. test_log_to_file() does this.
    """
    base_logger = logging.getLogger("astar")
    handlers = base_logger.handlers
    # Disable handlers
    base_logger.handlers = []
    # Make sure logging can reach pytest's caplog
    base_logger.propagate = True
    yield
    # Restore
    base_logger.handlers = handlers
    base_logger.propagate = False


@pytest.fixture(scope="package")
def mock_dir():
    """Path to mock directory."""
    return MOCK_DIR


@pytest.fixture(scope="package")
def mock_path():
    """Path to mock files."""
    return MOCK_DIR / "files"


@pytest.fixture(scope="class")
def patch_mock_path(mock_path):
    """Patch __search_path__ with test files mock path.

    Use only when needed internally, refer to filenames in tests using full
    absolute path (with the help of the mock_path fixture).
    """
    with patch("scopesim.rc.__search_path__", [mock_path]):
        yield


@pytest.fixture(scope="class")
def patch_all_mock_paths(mock_dir):
    with patch("scopesim.rc.__search_path__", UniqueList([mock_dir])):
        patched = {"!SIM.file.local_packages_path": str(mock_dir)}
        with patch.dict("scopesim.rc.__config__", patched):
            with patch.dict("scopesim.rc.__currsys__", patched):
                yield


@pytest.fixture(scope="package")
def mock_path_yamls():
    """Path to mock yaml files."""
    return MOCK_DIR / "yamls"


@pytest.fixture(scope="package")
def mock_path_micado():
    """Path to MICADO mock files."""
    return MOCK_DIR / "MICADO_SCAO_WIDE"


@pytest.fixture(scope="class")
def patch_mock_path_micado(mock_path_micado):
    """Patch __search_path__ with MICADO mock path.

    Use only when needed internally, refer to filenames in tests using full
    absolute path (with the help of the mock_path_micado fixture).
    """
    with patch("scopesim.rc.__search_path__", [mock_path_micado]):
        yield


@pytest.fixture(scope="function")
def no_file_error():
    """Patch currsys to avoid missing file error."""
    patched = {"!SIM.file.error_on_missing_file": False}
    with patch.dict("scopesim.rc.__currsys__", patched):
        yield


@pytest.fixture(scope="function")
def protect_currsys():
    """Prevent modification of global currsys."""
    with patch("scopesim.rc.__currsys__"):
        yield
