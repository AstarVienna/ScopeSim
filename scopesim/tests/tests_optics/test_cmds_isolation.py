# -*- coding: utf-8 -*-
"""Commands must not be shared between independently built objects.

These are the invariants that replaced the global ``rc.__currsys__``:
resolving a bang-string needs an explicit ``UserCommands``, anything built
without one gets its own rather than a process-wide dict, and ``cmds`` is
a named parameter everywhere so it can never land in ``.meta``.
"""

import pytest

from scopesim import rc
from scopesim.commands import UserCommands
from scopesim.effects import TERCurve
from scopesim.effects.effects import Effect
from scopesim.optics.optics_manager import OpticsManager
from scopesim.utils import from_currsys, default_cmds

# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring


class TestNoGlobalCurrsys:
    def test_rc_has_no_currsys(self):
        assert not hasattr(rc, "__currsys__")

    def test_from_currsys_requires_cmds(self):
        with pytest.raises(ValueError):
            from_currsys("!SIM.random.seed")

    def test_default_cmds_are_not_shared(self):
        one, two = default_cmds(), default_cmds()
        assert one is not two
        one["!TEL.area"] = 42
        assert two["!TEL.area"] != 42

    def test_writing_to_default_cmds_leaves_config_alone(self):
        before = rc.__config__["!TEL.area"]
        default_cmds()["!TEL.area"] = 99
        assert rc.__config__["!TEL.area"] == before


class TestEffectsGetOwnCmds:
    def test_effect_without_cmds_gets_its_own(self):
        one = TERCurve(wavelength=[1, 2], transmission=[1, 1])
        two = TERCurve(wavelength=[1, 2], transmission=[1, 1])
        assert one.cmds is not None
        assert one.cmds is not two.cmds

    def test_effects_do_not_leak_settings_to_each_other(self):
        one = TERCurve(wavelength=[1, 2], transmission=[1, 1])
        two = TERCurve(wavelength=[1, 2], transmission=[1, 1])
        one.cmds["!SIM.spectral.spectral_bin_width"] = 0.123
        assert from_currsys("!SIM.spectral.spectral_bin_width",
                            two.cmds) != 0.123

    def test_cmds_is_not_metadata(self):
        eff = TERCurve(wavelength=[1, 2], transmission=[1, 1],
                       cmds=UserCommands())
        assert "cmds" not in eff.meta


class TestSubclassGuard:
    """Effect.__init_subclass__ keeps cmds out of **kwargs for good."""

    def test_subclass_without_cmds_param_is_rejected(self):
        with pytest.raises(TypeError, match="explicit 'cmds' parameter"):
            class BadEffect(Effect):        # pylint: disable=unused-variable
                def __init__(self, **kwargs):
                    super().__init__(**kwargs)

    def test_subclass_with_cmds_param_is_accepted(self):
        class GoodEffect(Effect):
            def __init__(self, cmds=None, **kwargs):
                super().__init__(cmds=cmds, **kwargs)

        assert GoodEffect(cmds=UserCommands()).cmds is not None

    def test_subclass_without_own_init_is_accepted(self):
        class InheritingEffect(Effect):
            """Inherits an __init__ that already handles cmds."""

        assert InheritingEffect(cmds=UserCommands()).cmds is not None


class TestOpticsManagerDerivedParams:
    def test_requires_cmds(self):
        with pytest.raises(ValueError):
            OpticsManager()

    def test_derived_params_stay_out_of_the_shared_defaults(self):
        cmds = UserCommands()
        cmds["!INST.pixel_scale"] = 0.004
        cmds["!TEL.area"] = 100
        OpticsManager(cmds=cmds)
        # !TEL.etendue is computed and written into cmds ...
        assert cmds["!TEL.etendue"].value != 0
        # ... but must not reach the package defaults that every other
        # UserCommands chains onto. That leak is what rc.__currsys__ was.
        assert UserCommands()["!TEL.etendue"] == 0
        assert rc.__config__["!TEL.etendue"] == 0
