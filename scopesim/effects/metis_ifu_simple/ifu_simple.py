# -*- coding: utf-8 -*-
"""Effects for the METIS ifu_cube mode"""

from typing import ClassVar

from ..effects import Effect
from ...optics.fov import FieldOfView
from ...utils import from_currsys, get_logger

logger = get_logger(__name__)

class LineSpreadFunction(Effect):
    """
    Compute and apply line spread function to IFU cube
    """
    z_order: ClassVar[tuple[int, ...]] = (660,)
    report_plot_include: ClassVar[bool] = True
    report_table_include: ClassVar[bool] = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        params = {

        }
        self.meta.update(params)
        self.meta.update(kwargs)
        self.meta = from_currsys(self.meta, self.cmds)

        self.lsfwidth = None #self.get_width()
        self.kernel = None   #self.get_kernel()

    def apply_to(self, obj, **kwargs):
        """Apply the LSF"""
        if not isinstance(obj, FieldOfView):
            return obj

        lamc = from_currsys("!OBS.wavelen", self.cmds)
        print(">>>>>>>>> Central wavelength:", lamc)
