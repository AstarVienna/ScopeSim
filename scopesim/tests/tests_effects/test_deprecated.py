"""Test that calling deprecated effects results in errors"""
import pytest

from scopesim.effects.electronic import Quantization, SummedExposure

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

def test_error_on_quantization():
    with pytest.raises(AttributeError):
        Quantization()

def test_error_on_summed_exposure():
    with pytest.raises(AttributeError):
        SummedExposure()
