# -*- coding: utf-8 -*-
"""Contains much of the core functionality of ScopeSim."""

from .image_plane import ImagePlane
from . import image_plane_utils

from .optical_element import OpticalElement
from .optical_train import OpticalTrain
from .optics_manager import OpticsManager

from .surface import SpectralSurface
from . import surface_utils

from .fov import FieldOfView, FieldOfView1D, FieldOfView2D, FieldOfView3D
from .fov_manager import FOVManager
from .fov_volume_list import FovVolumeList
