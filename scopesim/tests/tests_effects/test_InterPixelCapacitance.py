# -*- coding: utf-8 -*-
"""Tests for InterPixelCapacitance effect."""

import pytest
import yaml
import numpy as np
from astropy.io import fits

from scopesim.effects.electronic import InterPixelCapacitance as IPC
from scopesim.detector import Detector

@pytest.fixture(scope="module", name="yaml_dict")
def fixture_yaml_dict():
    return yaml.full_load("""
    kernel: [[0., 0.02, 0.], [0.02, 0.92, 0.02], [0., 0.02, 0.]]
    """)


# pylint: disable=missing-class-docstring,missing-function-docstring
class TestInit:
    def test_initialises_correctly(self):
        ipc = IPC()
        assert isinstance(ipc, IPC)

    def test_initialised_with_nothing(self):
        ipc = IPC()
        assert np.all(ipc.kernel == np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]]))

    def test_initialises_with_kernel(self):
        kern = np.random.rand(3, 3)
        ipc = IPC(kernel=kern)
        assert np.all(ipc.kernel == kern)

    @pytest.mark.parametrize(
        "a_edge, a_corner, a_cross, kern",
        [(0, 0, 0, [[0, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]]),
         (0.02, 0, 0, [[0, 0.02, 0],
                       [0.02, 0.92, 0.02],
                       [0, 0.02, 0]]),
         (0.02, 0.002, 0, [[0.002, 0.02, 0.002],
                           [0.02, 0.912, 0.02],
                           [0.002, 0.02, 0.002]]),
         (0.0145, 0.0011, 0.0018, [[0.0011, 0.0127, 0.0011],
                                   [0.0163, 0.9376, 0.0163],
                                   [0.0011, 0.0127, 0.0011]])
         ])
    def test_initialises_with_params(self, a_edge, a_corner, a_cross, kern):
        ipc = IPC(alpha_edge=a_edge, alpha_corner=a_corner, alpha_cross=a_cross)
        assert np.allclose(ipc.kernel, np.asarray(kern))

    def test_initialises_from_yaml(self, yaml_dict):
        ipc = IPC(**yaml_dict)
        assert np.all(ipc.kernel == np.asarray(yaml_dict['kernel']))

    def test_refuse_negative_kernel(self):
        with pytest.raises(ValueError):
            _ = IPC(kernel=[[0, 0, 0], [0, -1, 0], [0, 0, 0]])

    def test_normalise_kernel_if_necessary(self):
        ipc = IPC(kernel = [[1., 1, 1], [1, 1, 1], [1, 1, 1]])
        assert np.allclose(ipc.kernel.sum(), 1)

    def test_str_shows_parameter_value(self):
        ipc = IPC(alpha_edge=0.02)
        assert "alpha_edge   = 0.02" in str(ipc)
        assert "alpha_corner = NA" in str(ipc)

@pytest.fixture(name="detector", scope="class")
def fixture_detector():
    """Instantiate a random detector object"""
    hdu = fits.ImageHDU(data=np.random.randn(214, 333))
    det = Detector(hdu.header)
    det._hdu.data = hdu.data
    return det

class TestApply:
    def test_reject_non_detector(self):
        ipc = IPC()
        obj = "Hallo"
        newobj = ipc.apply_to(obj)
        assert newobj == obj

    def test_convolve_unit_source(self):
        hdu = fits.ImageHDU(data=np.array([[0., 0., 0.],
                                           [0., 1., 0.],
                                           [0., 0., 0.]]))
        det = Detector(hdu.header)
        det._hdu.data = hdu.data
        ipc = IPC(kernel=np.random.rand(3, 3))
        newdet = ipc.apply_to(det)
        assert np.all(newdet.hdu.data == ipc.kernel)

    def test_preserves_shape(self, detector):
        oldshape = detector.hdu.data.shape
        ipc = IPC(kernel=np.random.rand(3, 3))
        newdet = ipc.apply_to(detector)
        newshape = newdet.hdu.data.shape
        assert newshape == oldshape

    def test_correlates_noise(self, detector):
        # IPC correlates pixel noise and reduces rms
        oldrms = np.std(detector.hdu.data - detector.hdu.data.mean())
        ipc = IPC(kernel=np.random.rand(3, 3))
        ipc.kernel /= np.sum(ipc.kernel)  # need to normalise for this test
        newdet = ipc.apply_to(detector)
        newrms = np.std(newdet.hdu.data - newdet.hdu.data.mean())
        assert newrms < oldrms
