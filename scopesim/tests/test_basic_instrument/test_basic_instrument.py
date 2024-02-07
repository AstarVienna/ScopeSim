import pytest
from pytest import raises
import os
from time import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u

import scopesim as sim
from scopesim.source import source_templates as st


PLOTS = False

inst_pkgs = os.path.join(os.path.dirname(__file__), "../mocks")
sim.rc.__currsys__["!SIM.file.local_packages_path"] = inst_pkgs


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
        spec = np.zeros(len(wave))
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=[0], y=[0], ref=[0], weight=[1e-3])

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["spectroscopy"])
        opt = sim.OpticalTrain(cmd)
        for effect_name in ["shot_noise", "dark_current", "readout_noise",
                            "atmospheric_radiometry", "source_fits_keywords",
                            "effects_fits_keywords", "config_fits_keywords"]:
            opt[effect_name].include = False

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
class TestObserveIfuMode:
    def test_runs(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.zeros(len(wave))
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        x = [-4, -2, 0, 2, 4] * 5
        y = [y0 for y0 in range(-4, 5, 2) for i in range(5)]
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=x, y=y, ref=[0]*len(x), weight=[1e-3]*len(x))

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["ifu"])
        opt = sim.OpticalTrain(cmd)
        for effect_name in ["shot_noise", "dark_current", "readout_noise",
                            "atmospheric_radiometry", "source_fits_keywords",
                            "effects_fits_keywords", "config_fits_keywords"]:
            opt[effect_name].include = False

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

        xs = [(157, 226), (317, 386), (476, 545), (640, 706), (797, 866)]
        spot_flux = 100000       # due to psf flux in large IFU slices (2")
        for i in range(5):
            x0, x1 = xs[i]
            trace_flux = det_im[:, x0:x1].sum()     # sum along a trace
            assert round(trace_flux / spot_flux) == 15 * 5

    def test_random_star_field(self):
        src = sim.source.source_templates.star_field(n=100, mmin=8, mmax=18, width=10)

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["ifu"])
        opt = sim.OpticalTrain(cmd)
        for effect_name in ["source_fits_keywords", "effects_fits_keywords",
                            "config_fits_keywords"]:
            opt[effect_name].include = False

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

    def test_runs_with_a_single_point_source(self):
        """
        Testing to see if point source halos can be seen in slit that is not
        centred on the object.

        Two ways of testing this by setting the ApertureList.meta values
        - Easy, computationally expensive :
          Make a FOV for each slit, then extend the FOV past the slit by x.
          This would be useful for MOS and LSS. E.g:
          - fov_for_each_aperture = True
          - extend_fov_beyond_slit = 2

        - Difficult, computationally cheap :
          Make a FOV that covers the whole IFU image slicer. E.g:
            - fov_for_each_aperture = False
            - extend_fov_beyond_slit = 0

        Current problems with each approach:
        - Easy : The whole FOV is added to the trace, not just the slit image
          Fix this by cutting a sub-fov from the fov for the XiLamImage
          [DONE]

        - Difficult : The FOVs are made by ApertureList THEN SpectralTraceList
          I'm not quite sure what SpTL is doing with the FovVolumes
          It calls the SpT.fov_grid method, which is a red flag.
          May need a re-write of SpT and SpTL

        """

        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.zeros(len(wave))
        spec[25::50] += 100      # every 0.05µm, offset by 0.025µm
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=[0], y=[0], ref=[0], weight=[1e-3])

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["ifu"])
        # Expand PSF to go over slit width (5x 2 arcsec -> IFU dims = 14" x 10")
        cmd["!OBS.psf_fwhm"] = 4
        # Create a single FOV for all 5 apertures -> False
        # cmd.yaml_dicts[3]["effects"][3]["kwargs"]["fov_for_each_aperture"] = True
        # FOV image extends past slit boundaries. Default = 0
        # cmd.yaml_dicts[3]["effects"][3]["kwargs"]["extend_fov_beyond_slit"] = 5

        opt = sim.OpticalTrain(cmd)
        for effect_name in ["shot_noise", "dark_current", "readout_noise",
                            "atmospheric_radiometry", "source_fits_keywords",
                            "effects_fits_keywords", "config_fits_keywords"]:
            opt[effect_name].include = False

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

        assert imp_im.sum() == pytest.approx(5251, rel=1e-3)


@pytest.mark.usefixtures("protect_currsys", "patch_all_mock_paths")
class TestObserveMosMode:
    def test_loads(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.ones(len(wave))
        spec[25::50] += 1      # every 0.05µm, offset by 0.025µm
        spec *= 1e3
        src = sim.Source(lam=wave*u.um, spectra=spec,
                         x=[-5, 5, 0, -5, 5],
                         y=[-5, -5, 0, 5, 5],
                         ref=[0]*5,
                         weight=[1e-3]*5)

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["mos"])
        opt = sim.OpticalTrain(cmd)
        opt.observe(src)
        hdul = opt.readout()[0]

        imp_im = opt.image_planes[0].data
        det_im = hdul[1].data

        if PLOTS:
            plt.subplot(121)
            plt.imshow(imp_im)
            plt.subplot(122)
            plt.imshow(det_im)
            plt.show()

    def test_photon_flux_remains_constant(self):
        wave = np.arange(0.7, 2.5, 0.001)
        spec = np.ones(len(wave))
        spec[25::50] += 1  # every 0.05µm, offset by 0.025µm
        spec *= 1e3
        src = sim.Source(lam=wave * u.um, spectra=spec,
                         x=[-5, 5, 0, -5, 5],
                         y=[-5, -5, 0, 5, 5],
                         ref=[0]*5,
                         weight=[1e-3]*5)

        cmd = sim.UserCommands(use_instrument="basic_instrument",
                               set_modes=["mos"])
        opt = sim.OpticalTrain(cmd)
        opt.observe(src)

        imp_im = opt.image_planes[0].data

        in_flux = 5*np.trapz(spec, wave)
        out_flux = np.sum(imp_im)

        if PLOTS:
            plt.imshow(imp_im)
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
        assert hdr["SIM EFF15 class"] == 'SourceDescriptionFitsKeywords'
        assert hdr["SIM CONFIG OBS filter_name"] == 'J'
        assert hdr["ESO ATM SEEING"] == opt.cmds["!OBS.psf_fwhm"]
