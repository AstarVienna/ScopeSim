# -*- coding: utf-8 -*-
"""Effects related to the conversion of photons into electrons.

   - LinearityCurve: Detector linearity
   - ADConversion:   Conversion from electrons to ADU

Related effects:
   - QuantumEfficiencyCurve: can be found in ter_curves.py
"""

from typing import ClassVar

import numpy as np

from .. import Effect
from ...detector import Detector
from ...utils import figure_factory, check_keys
from ...utils import from_currsys, real_colname
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
        if not isinstance(obj, Detector):
            return obj

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


class ADConversion(Effect):
    """Analogue-Digital Conversion effect.

    The effect applies the gain factor (electrons/ADU) to the detector readouts
    and converts the output to the desired data type (e.g. uint16).

    Examples
    --------
    The effect is usually instantiated in a yaml file.

    For a single-detector instrument:
    - name: ad_conversion
      description: Apply gain and convert electron count into integers
      class: ADConversion
      kwargs:
         dtype: uint16
         gain: "!DET.gain"       # or a number if !DET.gain has not been set yet

    For a multi-detector instrument the detector ids need to be identical to
    those used to instantiate the DetectorList effect.
    - name: ad_conversion
      description: Apply gain and convert electron count into integers
      class: ADConversion
      kwargs:
         dtype: uint16
         gain:
           id1:  2.2
           id2:  2.1
           id3   2.3

    Again, `!DET.gain` can be used here. This can be useful when the
    `DetectorModePropertiesSetter` effect is used to switch between different
    detector modes with different gain values.

    .. versionchanged:: 0.10.0

       Renamed from `Quantization` to `ADConversion`.
    """

    z_order: ClassVar[tuple[int, ...]] = (825,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {
            "dtype": "uint16",
            "gain": 1.      # default, usually overridden from yaml
        }
        self.meta.update(params)
        self.meta.update(kwargs)

    def _should_apply(self) -> bool:
        """Check cases where the effect should not be applied

        This does not do anything right now.
        """
        if self.cmds is None:
            logger.warning("Cannot access cmds for ADConversion effect.")
            return True

        # ..todo: need to deal with this case more realistically
        # Is this still necessary?
        #if self.cmds.get("!OBS.autoexpset", False):
        #    logger.info("DIT, NDIT determined by AutoExposure -> "
        #                "Create float32 output.")
        #    return False

        #if self.cmds["!OBS.ndit"] > 1:
        #    logger.info("NDIT larger than 1 -> Create float32 output.")
        #    return False

        return True

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, Detector):
            return obj

        new_dtype = self.meta["dtype"]

        # Apply the gain value (copy from DarkCurrent)
        # Note that this does not cater for the case where the gain is given
        # as a plain dictionary. Should we implement that?
        logger.info(f"Applying gain {from_currsys(self.cmds['!DET.gain'])}")
        if hasattr(self.cmds["!DET.gain"], "dic"):
            dtcr_id = obj.meta[real_colname("id", obj.meta)]
            gain = self.cmds["!DET.gain"].dic[dtcr_id]
        elif isinstance(self.cmds["!DET.gain"], (float, int)):
            gain = self.cmds["!DET.gain"]
        else:
            raise ValueError("<ADConversion>.meta['gain'] must be either "
                             f"dict or float, but is {self.cmds['!DET.gain']}")

        # Apply gain
        obj._hdu.data /= gain

        # Type-conversion wraps around input values that are higher or lower than
        # the respective maximum and minimum values of the new data type. Before
        # conversion to integer types with limited range (this is in particular
        # the case for 16 bits), we therefore need to cap the data.
        if np.issubdtype(new_dtype, np.integer):
            minval = np.iinfo(new_dtype).min
            maxval = np.iinfo(new_dtype).max
            minvals_mask = obj._hdu.data < minval
            maxvals_mask = obj._hdu.data > maxval
            if minvals_mask.any():
                obj._hdu.data[minvals_mask] = minval
                logger.warning(
                    f"Effect ADConversion: {minvals_mask.sum()} negative pixels")
            if maxvals_mask.any():
                obj._hdu.data[maxvals_mask] = maxval
                logger.warning(
                    f"Effect ADConversion: {maxvals_mask.sum()} saturated pixels")

        # This used to create a new ImageHDU with the same header but the data
        # set to the modified data. It should be fine to simply re-assign the
        # data attribute, but just in case it's not...
        logger.info("Applying digitization to dtype %s.", new_dtype)
        obj._hdu.data = obj._hdu.data.astype(new_dtype)

        return obj
