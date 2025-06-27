# -*- coding: utf-8 -*-
"""Effects for the METIS IFU_SMPL mode.

.. versionadded:: 0.10.0

Classes:
- LineSpreadFunction

"""

from .ifu_simple import LineSpreadFunction
from ...utils import get_logger

logger = get_logger(__name__)

__all__ = ["LineSpreadFunction",]
