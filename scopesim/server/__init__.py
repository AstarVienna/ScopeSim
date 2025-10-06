# -*- coding: utf-8 -*-
"""Utilities for remotely stored data.

Public functions
----------------
These functions are meant to be used by the end user.

- ``list_packages``: List available IRDB packages or versions of one package.
- ``download_packages``: Download and unpack IRDB package(s).
- ``list_example_data``: List available example datasets.
- ``download_example_data``: Download and cache example datasets.

Internal functions
------------------
These functions are meant for internal use elsewhere in the codebase, use only
if you know what you're doing and why you need them.

- ``create_retriever``: Create ``Pooch`` instance for handling cached downloads.
- ``check_packages``: Check if required package is in CWD, download if needed.
- ``get_all_packages_on_server``: Retrieve all unique package names present on
  server in known folders.

Exception classes
-----------------

- ``ServerError``: Error in accessing the server, usually raised from various
  other errors further down.
- ``PkgNotFoundError``: Specified IRDB package or specific release of an IRDB
  package wasn't found on the server.

"""

from .download_utils import ServerError, create_retriever
from .database import (
    PkgNotFoundError,
    download_packages,
    list_packages,
    get_all_packages_on_server,
    check_packages,
)
from .example_data_utils import download_example_data, list_example_data
