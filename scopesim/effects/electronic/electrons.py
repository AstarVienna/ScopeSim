# -*- coding: utf-8 -*-
"""Effects related to the conversion of photons into electrons.

   - LinearityCurve: Detector linearity
   - ADConversion:   Conversion from electrons to ADU

Related effects:
   - QuantumEfficiencyCurve: can be found in ter_curves.py
"""

from typing import ClassVar

import numpy as np
from scipy.signal import oaconvolve

from .. import Effect
from ...detector import Detector
from ...utils import figure_factory, check_keys
from ...utils import from_currsys, real_colname, get_logger
from . import logger

logger = get_logger(__name__)

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


class InterPixelCapacitance(Effect):
    r"""Inter-pixel capacitance effect.

    The effect models cross-talk due to inter-pixel capacitance with
    a convolution kernel following [1]_.

    .. versionadded:: 0.11.1

    Example
    -------
    The effect is usually instantiated in a yaml file.

    The first example uses the three-parameter model in Eq. (9) of [1].
    The three parameters are `alpha_edge` (corresponding to :math:`\alpha`) for
    the effect of neighbouring pixels sharing an edge with the pixel under
    consideration, `alpha_corner` (corresponding to :math:`\alpha^\prime`) for
    pixels sharing a corner, and `alpha_aniso` (corresponding to
    :math:`\alpha_{+}`) to allow for different capacitive coupling along rows
    and columns. The simpler one- and two-parameters models are recovered by
    setting `alpha_aniso` and/or `alpha_corner` to zero.

    ::

      - name: ipc
        description: Apply inter-pixel capacitance
        class: InterPixelCapacitance
        kwargs:
           alpha_edge: 0.02
           alpha_corner: 0.002
           alpha_aniso: 0

    Alternatively, a convolution kernel can be provided explicitely:

    ::

      - name: ipc
        description: Apply inter-pixel capacitance
        class: InterPixelCapacitance
        kwargs:
           kernel: [
              [0.0011, 0.0127, 0.0011],
              [0.0163, 0.9360, 0.0164],
              [0.0011, 0.0127, 0.0011],
            ]

    .. [1] Kannawadi et al., "The Impact of Interpixel Capacitance in CMOS
       Detectors on PSF Shapes and Implications for WFIRST",
       PASP 128, 095001 (2016).
    """

    z_order: ClassVar[tuple[int, ...]] = (810,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.meta.update(kwargs)
        self.kernel = self._build_kernel(kwargs)

    def _build_kernel(self, params):
        """Build a 3x3 kernel."""
        if "kernel" in params:
            kernel = np.asarray(params['kernel']).astype(float)
            kernsum = np.sum(kernel)
            if kernsum > 1:
                logger.warning("IPC kernel is larger than one, normalising")
                kernel /= kernsum
            if kernsum <= 0:
                raise ValueError("IPC kernel has negative normalisation")
            return kernel

        a_corner = params.get("alpha_corner", 0)
        a_aniso = params.get("alpha_aniso", 0)
        a_edge = params.get("alpha_edge", 0)

        kernel = np.array([
            [a_corner, a_edge - a_aniso, a_corner],
            [a_edge + a_aniso, 1 - 4 * (a_edge + a_corner), a_edge + a_aniso],
            [a_corner, a_edge - a_aniso, a_corner]
        ])
        return kernel

    def apply_to(self, det, **kwargs):
        if not isinstance(det, Detector):
            logger.debug("%s applied to %s", self.display_name,
                         det.__class__.__name__)
            return det

        newdata = oaconvolve(det._hdu.data, self.kernel, mode="same")
        det._hdu.data = newdata
        return det

    def update(self, **kwargs):
        """Update the IPC kernel.

        An instance `ipc` of `InterPixelCapacitance` can be updated by
        specifying either a new `kernel`:
        ``ipc.update(kernel=[[0., 0.02, 0], [0, 0.92, 0], [0, 0.02, 0]])``
        or by specifying one or more of the `alpha` parameters:
        ``ipc.update(alpha_edge=0.02, alpha_corner=0.002)``

        Notes
        -----
        Unspecified `alpha` parameters (here `alpha_aniso`) default to zero.
        """
        if "kernel" in kwargs:
            for key in ["alpha_edge", "alpha_corner", "alpha_aniso"]:
                self.meta.pop(key, None)
        self.meta.update(kwargs)
        self.kernel = self._build_kernel(kwargs)

    def __str__(self):
        kernel_str = np.array2string(
            self.kernel,
            precision=4,
            floatmode="fixed",
            prefix="   kernel = ",
        )
        msg = (
            f"<{self.__class__.__name__}> \"{self.meta['description']}\" :"
            f"alpha_edge   = {self.meta.get('alpha_edge', 'NA')}"
            f"alpha_corner = {self.meta.get('alpha_corner', 'NA')}"
            f"alpha_aniso  = {self.meta.get('alpha_aniso', 'NA')}"
            f"kernel = {kernel_str}"
        )
        return msg


class ADConversion(Effect):
    """Analogue-Digital Conversion effect.

    The effect applies the gain factor (electrons/ADU) to the detector readouts
    and converts the output to the desired data type (e.g. uint16).

    Examples
    --------
    The effect is usually instantiated in a yaml file.

    For a single-detector instrument:

    ::

      - name: ad_conversion
        description: Apply gain and convert electron count into integers
        class: ADConversion
        kwargs:
           dtype: uint16
           gain: "!DET.gain"  # or a number if !DET.gain has not been set yet

    For a multi-detector instrument the detector ids need to be identical to
    those used to instantiate the DetectorList effect.

    ::

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
        """Check cases where the effect should not be applied.

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
