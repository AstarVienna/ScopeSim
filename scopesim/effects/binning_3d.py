# -*- coding: utf-8 -*-
"""Where else to put this?."""

from typing import ClassVar

from astropy import units as u

from .effects import Effect
from ..optics.fov import FieldOfView
from ..utils import get_logger, from_currsys, unit_includes_per_physical_type


logger = get_logger(__name__)


class FluxBinning3D(Effect):
    """Takes care of cube flux conversion in absence of a SpectralTraceList.

    .. versionadded:: 0.10.0

    """

    z_order: ClassVar[tuple[int, ...]] = (690,)

    def apply_to(self, fov, **kwargs):
        """See parent docstring."""
        if not isinstance(fov, FieldOfView):
            return fov

        if fov.hdu is None or fov.hdu.header["NAXIS"] != 3:
            logger.error("Cannot apply cube flux binning.")
            return fov

        bunit = u.Unit(fov.hdu.header["BUNIT"].lower())
        pixel_area = fov.pixel_area << u.arcsec**2

        # Spatial binning
        if unit_includes_per_physical_type(bunit, "solid angle"):
            fov.hdu.data *= pixel_area.value
            bunit *= pixel_area.unit
        else:
            logger.warning("Cube is already binned spatially.")

        # Spectral binning
        if unit_includes_per_physical_type(bunit, "length"):
            fov.hdu.data *= self.dwave.value
            bunit *= self.dwave.unit
        else:
            logger.warning("Cube is already binned spectrally.")

        fov.hdu.header["BUNIT"] = bunit.to_string("fits")

        # This is done in SpectralTraceList as well, idk...
        fov.cube = fov.hdu
        return fov

    @property
    def dwave(self) -> u.Quantity[u.um]:
        # TODO: turn into class attribute once cmds lookup works...
        dwave = from_currsys("!SIM.spectral.spectral_bin_width", self.cmds)
        return dwave << u.um
