# -*- coding: utf-8 -*-
"""Effects related to the conversion of photons into electrons.

   - LinearityCurve: Detector linearity
   - Quantization:   Conversion from electrons to ADU

Related effects:
   - QuantumEfficiencyCurve: can be found in ter_curves.py
"""

from typing import ClassVar

import numpy as np

from .. import Effect
from ...base_classes import DetectorBase
from ...utils import from_currsys, figure_factory, check_keys
from . import logger


class LinearityCurve(Effect):
    """
    Detector linearity effect.

    The detector linearity curve is set in terms of `incident` flux (e/s) and
    `measured` detector values (ADU).

    Examples
    --------
    The effect can be instantiated in various ways.::

        - name: detector_linearity
          class: LinearityCurve
          kwargs:
            filename: FPA_linearity.dat

        - name: detector_linearity
          class: LinearityCurve
          kwargs:
            array_dict: {incident: [0, 77000, 999999999999],
                         measured: [0, 77000, 77000]}

        - name: detector_linearity
          class: LinearityCurve
          kwargs:
            incident: [0, 77000, 99999999]
            measured: [0, 77000, 77000]

    """

    required_keys = {"ndit"}
    z_order: ClassVar[tuple[int, ...]] = (840,)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.meta.update(kwargs)

        check_keys(self.meta, self.required_keys, action="error")

    def apply_to(self, obj, **kwargs):
        if isinstance(obj, DetectorBase):
            ndit = from_currsys(self.meta["ndit"], self.cmds)
            if self.table is not None:
                incident = self.table["incident"] * ndit
                measured = self.table["measured"] * ndit
            else:
                incident = np.asarray(from_currsys(self.meta["incident"],
                                                   self.cmds)) * ndit
                measured = np.asarray(from_currsys(self.meta["measured"],
                                                   self.cmds)) * ndit
            obj._hdu.data = np.interp(obj._hdu.data, incident, measured)

        return obj

    def plot(self, **kwargs):
        fig, ax = figure_factory()

        ndit = from_currsys(self.meta["ndit"], self.cmds)
        incident = self.table["incident"] * ndit
        measured = self.table["measured"] * ndit

        ax.loglog(incident, measured, **kwargs)
        ax.set_xlabel("Incident [ph s$^-1$]")
        ax.set_ylabel("Measured [e- s$^-1$]")

        return fig


class Quantization(Effect):
    """Converts raw data to whole photons."""

    z_order: ClassVar[tuple[int, ...]] = (825,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "dtype": "uint32",
        }
        self.meta.update(params)
        self.meta.update(kwargs)

    def _should_apply(self) -> bool:
        if self.cmds is None:
            logger.warning("Cannot access cmds for Quantization effect.")
            return True

        if self.cmds.get("!OBS.autoexpset", False):
            logger.debug("DIT, NDIT determined by AutoExposure -> "
                         "quantization is not applied.")
            return False

        if self.cmds["!OBS.ndit"] > 1:
            logger.debug("NDIT set to 1 -> quantization is not applied.")
            return False

        return True

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, DetectorBase):
            return obj

        if not self._should_apply():
            return obj

        new_dtype = self.meta["dtype"]
        if not np.issubdtype(new_dtype, np.integer):
            logger.warning("Setting quantized data to dtype %s, which is not "
                           "an integer subtype.", new_dtype)

        # Remove values below 0, because for unsigned types these get wrapped
        # around to MAXINT, which is even worse than negative values.
        # TODO: Write a test for this.
        if np.issubdtype(new_dtype, np.unsignedinteger):
            negvals_mask = obj._hdu.data < 0
            if negvals_mask.any():
                logger.warning(f"Effect Quantization: {negvals_mask.sum()} negative pixels")
                obj._hdu.data[negvals_mask] = 0

        # This used to create a new ImageHDU with the same header but the data
        # set to the modified data. It should be fine to simply re-assign the
        # data attribute, but just in case it's not...
        logger.debug("Applying quantization to dtype %s.", new_dtype)
        obj._hdu.data = np.floor(obj._hdu.data).astype(new_dtype)

        return obj
