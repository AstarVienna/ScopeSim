# -*- coding: utf-8 -*-
import pytest

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u

import scopesim as sim
from scopesim.source import source_templates as st
from scopesim.effects import AutoExposure, Quantization


PLOTS = False

SWITCHOFF = [
    "shot_noise",
    "dark_current",
    "readout_noise",
    "atmospheric_radiometry",
    "source_fits_keywords",
    "effects_fits_keywords",
    "config_fits_keywords",
]


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestLoadsUserCommands:
    def test_loads(self):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        assert isinstance(cmd, sim.UserCommands)
        assert cmd["!INST.pixel_scale"] == 0.2


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestLoadsOpticalTrain:
    def test_loads(self):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        opt = sim.OpticalTrain(cmd)
        assert isinstance(opt, sim.OpticalTrain)
        assert opt["#slit_wheel.current_slit"] == "!OBS.slit_name"
        assert opt["#slit_wheel.current_slit!"] == "narrow"


@pytest.mark.slow
@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestObserveImagingMode:
    def test_runs(self):
        src = st.star(flux=9)
        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["imaging"])
        opt = sim.OpticalTrain(cmd)
        opt.observe(src)
        hdul = opt.readout()[0]

        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(opt.image_planes[0].data, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im, norm=LogNorm())
            plt.show()

        assert det_im.max() > 1e5
        assert det_im[505:520, 505:520].sum() > 3e6


@pytest.mark.slow
@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestObserveSpectroscopyMode:
    """
    Test the number of spots along the three spectral traces.

    Spots are places at 0.05um intervals, offset by 0.025um

    3 traces vertically covering the regions:
    - J: 1.0, 1.35     dwave = 0.35 --> R~2800      -> 7 spots
    - H: 1.3, 1.8      dwave = 0.5  --> R~2000      -> 10 spots
    - K: 1.75, 2.5     dwave = 0.75 --> R~1300      -> 15 spots
    """

    def test_runs(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.zeros_like(wave)
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=[0], y=[0], ref=[0], weight=[1e-3])

        cmd = sim.UserCommands(
            use_instrument="basic_instrument",
            set_modes=["spectroscopy"],
            ignore_effects=SWITCHOFF,
        )
        opt = sim.OpticalTrain(cmd)

        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im)
            plt.show()

        xs = [slice(172, 199), slice(495, 525), slice(821, 850)]
        dlams = np.array([0.35, 0.5, 0.75])
        n_spots = (dlams / 0.05).round().astype(int)

        spot_flux = 28000       # due to psf flux losses in narrow slit (0.5")
        for sl, n in zip(xs, n_spots):
            trace_flux = det_im[:, sl].sum()     # sum along a trace
            assert round(trace_flux / spot_flux) == n


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestSourceImageNotAffected:
    """
    Test that an ImageHDU source object is not altered during a (spectroscopic) observation
    where the source is applied to multiple FOVs via FieldOfView._make_cube_imagefields()
    """

    @pytest.mark.slow
    def test_runs(self):
        src = st.uniform_source()
        src_int = np.mean(src.fields[0].field.data) # original source image intensity

        cmd = sim.UserCommands(
            use_instrument="basic_instrument",
            set_modes=["spectroscopy"],
            ignore_effects=SWITCHOFF,
        )
        opt = sim.OpticalTrain(cmd)

        opt.observe(src)

        # iterated through multiple FOVs?
        assert len(opt.fov_manager.fovs) > 1
        # source and fov.field objects not altered by observe()?
        assert np.mean(src.fields[0].field.data) == src_int
        assert np.mean(opt.fov_manager.fovs[0].fields[0].data) == src_int


@pytest.mark.slow
@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestObserveIfuMode:
    def test_runs(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.zeros_like(wave)
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        x = np.tile(np.arange(-4, 5, 2), 5)
        y = np.repeat(np.arange(-4, 5, 2), 5)
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=x, y=y, ref=[0]*len(x), weight=[1e-3]*len(x))

        cmd = sim.UserCommands(
            use_instrument="basic_instrument",
            set_modes=["ifu"],
            ignore_effects=SWITCHOFF,
        )
        opt = sim.OpticalTrain(cmd)

        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im)
            plt.show()

        xs = [slice(157, 226), slice(317, 386), slice(476, 545),
              slice(640, 706), slice(797, 866)]
        spot_flux = 100000       # due to psf flux in large IFU slices (2")
        for sl in xs:
            trace_flux = det_im[:, sl].sum()     # sum along a trace
            assert round(trace_flux / spot_flux) == 15 * 5

    def test_random_star_field(self):
        src = sim.source.source_templates.star_field(
            n=100, mmin=8, mmax=18, width=10)

        cmd = sim.UserCommands(
            use_instrument="basic_instrument",
            set_modes=["ifu"],
            ignore_effects=[
                "source_fits_keywords",
                "effects_fits_keywords",
                "config_fits_keywords",
            ],
        )
        opt = sim.OpticalTrain(cmd)

        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im, norm=LogNorm())
            plt.subplot(122)
            plt.imshow(det_im, norm=LogNorm())
            plt.show()


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestFitsHeader:
    def test_source_keywords_in_header(self):
        src = st.star()
        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["imaging"])
        opt = sim.OpticalTrain(cmd)
        opt.observe(src)
        hdul = opt.readout()[0]
        hdr = hdul[0].header

        assert hdr["SIM SRC0 object"] == 'star'
        assert hdr["SIM EFF14 class"] == 'SourceDescriptionFitsKeywords'
        assert hdr["SIM CONFIG OBS filter_name"] == 'J'
        assert hdr["ESO ATM SEEING"] == opt.cmds["!OBS.psf_fwhm"]


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestModeStatus:
    def test_concept_mode_init(self):
        with pytest.raises(NotImplementedError):
            _ = sim.UserCommands(use_instrument="basic_instrument",
                                 set_modes=["mock_concept_mode"])

    def test_concept_mode_change(self):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        with pytest.raises(NotImplementedError):
            cmd.set_modes("mock_concept_mode")

    def test_experimental_mode_init(self, caplog):
        _ = sim.UserCommands(use_instrument="basic_instrument",
                             set_modes=["mock_experimental_mode"])
        assert ("Mode 'mock_experimental_mode' is still in experimental stage"
                in caplog.text)

    def test_experimental_mode_change(self, caplog):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        cmd.set_modes("mock_experimental_mode")
        assert ("Mode 'mock_experimental_mode' is still in experimental stage"
                in caplog.text)

    def test_deprecated_mode_init(self):
        with pytest.raises(
                DeprecationWarning,
                match="Instrument mode 'mock_deprecated_mode' is deprecated."):
            _ = sim.UserCommands(use_instrument="basic_instrument",
                                 set_modes=["mock_deprecated_mode"])

    def test_deprecated_mode_change(self):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        with pytest.raises(
                DeprecationWarning,
                match="Instrument mode 'mock_deprecated_mode' is deprecated."):
            cmd.set_modes("mock_deprecated_mode")

    def test_deprecated_msg_mode_init(self):
        with pytest.raises(
                DeprecationWarning, match="This mode is deprecated."):
            _ = sim.UserCommands(use_instrument="basic_instrument",
                                 set_modes=["mock_deprecated_mode_msg"])

    def test_deprecated_msg_mode_change(self, caplog):
        cmd = sim.UserCommands(use_instrument="basic_instrument")
        with pytest.raises(
                DeprecationWarning, match="This mode is deprecated."):
            cmd.set_modes("mock_deprecated_mode_msg")


@pytest.fixture(scope="function", name="obs")
def basic_opt_observed():
    src = st.star(flux=15)
    cmd = sim.UserCommands(use_instrument="basic_instrument",
                           ignore_effects=SWITCHOFF,
                           properties={"!OBS.dit": 10, "!OBS.ndit": 1})
    opt = sim.OpticalTrain(cmd)
    opt.observe(src)
    default = int(opt.readout()[0][1].data.sum())

    quanteff = Quantization(cmds=opt.cmds)
    opt.optics_manager["basic_detector"].effects.append(quanteff)
    return opt, default, quanteff


@pytest.fixture(scope="function", name="obs_aeq")
def basic_opt_with_autoexp_and_quant_observed():
    src = st.star(flux=15)
    cmd = sim.UserCommands(use_instrument="basic_instrument",
                           ignore_effects=SWITCHOFF,
                           properties={"!OBS.dit": 10, "!OBS.ndit": 1})
    opt = sim.OpticalTrain(cmd)
    opt.observe(src)
    default = int(opt.readout()[0][1].data.sum())

    autoexp = AutoExposure(cmds=opt.cmds, mindit=1,
                           full_well=100, fill_frac=0.8)
    quanteff = Quantization(cmds=opt.cmds)
    opt.optics_manager["basic_detector"].effects.insert(2, autoexp)
    opt.optics_manager["basic_detector"].effects.append(quanteff)
    opt.cmds["!OBS.exptime"] = 60
    return opt, default, quanteff


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestDitNdit:
    @pytest.mark.parametrize(("dit", "ndit", "factor", "quant"),
                             [(20, 1, 2, True), (10, 3, 3, False)])
    def test_obs_dict(self, obs, dit, ndit, factor, quant):
        """This should just use dit, ndit from !OBS."""
        opt, default, quanteff = obs
        # This should work with patch.dict, but doesn't :(
        o_dit, o_ndit = opt.cmds["!OBS.dit"], opt.cmds["!OBS.ndit"]
        opt.cmds["!OBS.dit"] = dit
        opt.cmds["!OBS.ndit"] = ndit
        kwarged = int(opt.readout(reset=False)[0][1].data.sum())
        assert quanteff._should_apply() == quant
        opt.cmds["!OBS.dit"] = o_dit
        opt.cmds["!OBS.ndit"] = o_ndit
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    @pytest.mark.parametrize(("dit", "ndit", "factor", "quant"),
                             [(20, 1, 2, True),
                              (10, 3, 3, False),
                              (None, None, 1, True)])
    def test_kwargs_override_obs_dict(self, obs, dit, ndit, factor, quant):
        """This should prioritize kwargs and fallback to !OBS."""
        opt, default, quanteff = obs
        kwarged = int(opt.readout(dit=dit, ndit=ndit)[0][1].data.sum())
        assert quanteff._should_apply() == quant
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    @pytest.mark.parametrize(("dit", "ndit", "factor", "quant"),
                             [(20, 1, 2, True),
                              (10, 3, 3, False),
                              (None, None, 1, True)])
    def test_kwargs_override_obs_dict_also_with_autoexp(
            self, obs_aeq, dit, ndit, factor, quant):
        """This should prioritize dit, ndit from kwargs.

        Lacking those, dit and ndit from !OBS should be used over exptime.
        """
        opt, default, quanteff = obs_aeq
        kwarged = int(opt.readout(dit=dit, ndit=ndit, reset=False)[0][1].data.sum())
        assert quanteff._should_apply() == quant
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    @pytest.mark.parametrize(("exptime", "factor"),
                             [(20, 2), (30, 3), (None, 6)])
    def test_autoexp(self, obs_aeq, exptime, factor):
        """This should prioritize kwargs and fallback to !OBS."""
        opt, default, quanteff = obs_aeq
        # This should work with patch.dict, but doesn't :(
        o_dit, o_ndit = opt.cmds["!OBS.dit"], opt.cmds["!OBS.ndit"]
        opt.cmds["!OBS.dit"] = None
        opt.cmds["!OBS.ndit"] = None
        kwarged = int(opt.readout(exptime=exptime, reset=False)[0][1].data.sum())
        assert not quanteff._should_apply()
        opt.cmds["!OBS.dit"] = o_dit
        opt.cmds["!OBS.ndit"] = o_ndit
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    @pytest.mark.parametrize(("exptime", "factor", "quant"),
                             [(30, 3, False), (None, 1, True)])
    def test_autoexp_overrides_obs_dict(self, obs_aeq, exptime, factor, quant):
        """This should prioritize kwargs and use dit, ndit when None."""
        opt, default, quanteff = obs_aeq
        kwarged = int(opt.readout(exptime=exptime, reset=False)[0][1].data.sum())
        assert quanteff._should_apply() == quant
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    @pytest.mark.parametrize(("dit", "ndit", "factor", "quant"),
                             [(90, 1, 9, True), (2, 90, 18, False)])
    def test_ditndit_in_kwargs_while_also_having_autoexp(
            self, obs_aeq, dit, ndit, factor, quant):
        """This should prioritize dit, ndit from kwargs."""
        opt, default, quanteff = obs_aeq
        kwarged = int(opt.readout(dit=dit, ndit=ndit, reset=False)[0][1].data.sum())
        assert quanteff._should_apply() == quant
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    @pytest.mark.parametrize(("dit", "ndit", "exptime", "factor", "quant"),
                             [(90, 1, None, 9, True),
                              (2, 90, 20, 18, False)])
    def test_ditndit_in_kwargs_while_also_having_autoexp_and_exptime(
            self, obs_aeq, dit, ndit, exptime, factor, quant):
        """This should prioritize dit, ndit from kwargs and ignore exptime."""
        opt, default, quanteff = obs_aeq
        kwarged = int(opt.readout(exptime=exptime,
                                  dit=dit, ndit=ndit,
                                  reset=False)[0][1].data.sum())
        assert quanteff._should_apply() == quant
        # Quantization results in ~4% loss, which is fine:
        assert pytest.approx(kwarged / default, rel=.05) == factor

    def test_throws_for_no_anything(self, obs):
        """No specification whatsoever, so throw error."""
        opt, default, quanteff = obs
        opt.cmds["!OBS.exptime"] = None
        # This should work with patch.dict, but doesn't :(
        o_dit, o_ndit = opt.cmds["!OBS.dit"], opt.cmds["!OBS.ndit"]
        opt.cmds["!OBS.dit"] = None
        opt.cmds["!OBS.ndit"] = None
        with pytest.raises(ValueError):
            opt.readout(reset=False)
        opt.cmds["!OBS.dit"] = o_dit
        opt.cmds["!OBS.ndit"] = o_ndit

    def test_throws_for_no_ditndit_no_autoexp_kwargs(self, obs):
        """This should use exptime from kwargs, but fail w/o AutoExp."""
        opt, default, quanteff = obs
        opt.cmds["!OBS.exptime"] = None
        with pytest.raises(ValueError):
            opt.readout(exptime=60, reset=False)

    def test_throws_for_no_ditndit_no_autoexp_obs(self, obs):
        """This should fallback to !OBS.exptime, but fail w/o AutoExp."""
        opt, default, quanteff = obs
        opt.cmds["!OBS.exptime"] = 60
        # This should work with patch.dict, but doesn't :(
        o_dit, o_ndit = opt.cmds["!OBS.dit"], opt.cmds["!OBS.ndit"]
        opt.cmds["!OBS.dit"] = None
        opt.cmds["!OBS.ndit"] = None
        with pytest.raises(ValueError):
            opt.readout(reset=False)
        opt.cmds["!OBS.dit"] = o_dit
        opt.cmds["!OBS.ndit"] = o_ndit
