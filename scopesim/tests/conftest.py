# -*- coding: utf-8 -*-
"""Global fixtures for pytest."""

from pathlib import Path

import pytest
from unittest.mock import patch

import scopesim as sim

MOCK_DIR = Path(__file__).parent / "mocks"

sim.rc.__currsys__["!SIM.file.error_on_missing_file"] = True


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
    # FIXME: remove the second part of this asap
    try:
        with patch.dict("scopesim.rc.__currsys__", patched):
            yield
    except KeyError:
        with patch.dict("scopesim.rc.__currsys__.cmds", patched):
            yield