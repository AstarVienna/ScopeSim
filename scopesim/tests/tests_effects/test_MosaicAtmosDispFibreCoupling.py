import pytest
from pytest import approx

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u

from scopesim.effects import MosaicAtmosDispFibreCoupling
from scopesim import rc

PLOTS = False


class TestInit:
    def test_initialises_with_nothing(self):
        fc = MosaicAtmosDispFibreCoupling()
        assert isinstance(fc, MosaicAtmosDispFibreCoupling)

    def test_get_output_from_fibre_coupling_class(self):
        """
        Not sure if this is a worthy test, but it shows me what I want to see
        """

        fc = MosaicAtmosDispFibreCoupling()
        wave, trans = fc.coupling_fractions.run(0, 0, 0, "LR_NIR_H", 0.01, -1, -1)
        for tran in trans:
            assert np.all(tran > 0)


class TestApplyTo:
    pass
