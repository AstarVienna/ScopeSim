ScopeSim concepts and architecture
====================================

This page explains the key objects in ScopeSim and how they fit together.
Reading it will help you understand what happens during a simulation and how
to customise it.

---

The simulation pipeline
------------------------

Every ScopeSim simulation follows the same pipeline:

.. code-block:: text

    Source  ──►  OpticalTrain  ──►  Detector  ──►  FITS output
       ↑               ↑
   (target)     (telescope + instrument +
                   atmosphere + sky)

1. A **Source** object describes the on-sky target: spatial positions and
   spectra.
2. The **OpticalTrain** applies a sequence of **Effects** — PSFs, transmission
   curves, noise sources, detector geometry — that transform the source signal
   as photons travel from the sky to the detector.
3. The **Detector** integrates the focal-plane image, adds readout noise and
   dark current, and produces the final pixel array.
4. The result is returned as an ``astropy.fits.HDUList``.

The ``Simulation`` class is a convenience wrapper that bundles steps 2–4
together and reads the instrument configuration from the IRDB packages.

---

The five core objects
----------------------

Source
~~~~~~

A ``Source`` represents an on-sky target. It holds:

- **Spatial information** — positions of emitting elements (point sources or
  images)
- **Spectral information** — one or more spectra associated with those elements
- **Flux scaling** — brightnesses in physical units or magnitudes

Sources are created using ``scopesim_templates`` functions for common
astronomical objects, or directly from FITS images, tables, or spectra::

    import scopesim_templates as st

    # Point sources
    src = st.stellar.star(mag=20, filter_name="K")
    src = st.stellar.cluster(mass=1e4, distance=50000)

    # Extended source from a FITS image
    from scopesim import Source
    src = Source(image_hdu=my_fits_image_hdu, spectra=my_spectrum)

    # Combine sources with +
    combined = star_src + galaxy_src

OpticalTrain
~~~~~~~~~~~~~

The ``OpticalTrain`` is the heart of ScopeSim. It models everything between
the sky and the detector: atmosphere, telescope mirrors, AO system, instrument
optics, filters, PSF, and detectors. It contains an ordered list of
**Effects** that are applied in sequence when ``observe()`` is called.

You rarely create an ``OpticalTrain`` directly — the ``Simulation`` class does
it for you. But you can inspect and customise it::

    sim = Simulation("MICADO", mode=["MCAO_4mas", "IMG"])
    ot = sim.optical_train

    ot.effects                          # view all effects
    ot["detector_linearity"].include = False  # disable an effect
    ot.image_planes[0].data             # access intermediate focal-plane image

Effect
~~~~~~~

Effects are the building blocks of the optical model. Each ``Effect`` subclass
models one physical process:

- **SurfaceList** — mirrors and optical surfaces (transmission, emission)
- **FieldConstantPSF**, **FieldVaryingPSF** — PSF from a FITS cube
- **AnisocadoConstPSF** — AnisoCADO-based AO PSF
- **SeeingPSF**, **GaussianDiffractionPSF** — analytical PSFs
- **FilterCurve**, **QuantumEfficiencyCurve** — spectral transmission/efficiency
- **ShotNoise**, **DarkCurrent**, **ReadoutNoise** — detector noise sources
- **DetectorList** — detector geometry, pixel scale, gaps
- **AutoExposure**, **SummedExposure** — exposure logic

Effects are defined in the instrument YAML files in the IRDB. You can also
write custom effects — see the :doc:`examples/3_custom_effects` notebook.

UserCommands / settings
~~~~~~~~~~~~~~~~~~~~~~~~

The ``UserCommands`` object (accessible as ``sim.settings``) is a nested
dictionary that holds all simulation parameters. It is constructed from the
instrument YAML files but any parameter can be overridden at runtime.

Parameters are accessed using **bang-string** (``!``) notation::

    sim.settings["!OBS.dit"] = 60           # integration time [s]
    sim.settings["!OBS.ndit"] = 10          # number of integrations
    sim.settings["!OBS.filter_name"] = "K"  # filter selection
    sim.settings["!ATMO.seeing"] = 0.7      # seeing FWHM [arcsec]
    sim.settings["!SIM.spectral.wave_min"] = 1.9  # wavelength range [µm]

The ``!`` prefix and dot notation resolve paths through nested YAML
dictionaries. The top-level aliases are:

.. list-table::
   :widths: 15 85
   :header-rows: 1

   * - Alias
     - Covers
   * - ``!OBS``
     - Observation parameters (DIT, NDIT, filter, mode)
   * - ``!SIM``
     - Simulation parameters (wavelength range, pixel oversampling, file paths)
   * - ``!ATMO``
     - Atmospheric parameters (seeing, background, transmission)
   * - ``!TEL``
     - Telescope parameters (collecting area, emissivity)
   * - ``!INST``
     - Instrument parameters (pixel scale, plate scale)
   * - ``!DET``
     - Detector parameters (readout mode, gain, dark current)

Detector
~~~~~~~~~

The ``Detector`` object manages the focal-plane array. It holds one or more
``DetectorWindow`` instances (individual chips), integrates signal from the
focal plane, and applies per-chip noise models and non-linearity corrections.
You typically interact with the output HDUList rather than the ``Detector``
directly.

---

Fields of View (FOVs)
----------------------

When ``OpticalTrain.observe()`` is called, ScopeSim divides the simulation
volume into small *Fields of View* (FOVs). Each FOV covers a spatial
sub-region and a wavelength slice, chosen to match the spectral resolution of
the current setup. Effects are applied independently within each FOV, then
the results are assembled into the full focal-plane image.

This design enables:

- Efficient memory usage for large detector arrays
- Wavelength-dependent PSFs (field-varying PSF cubes)
- Multi-chip detectors with gaps and offsets
- Spectroscopic and IFU modes with curved/tilted traces

You rarely need to interact with FOVs directly. If you need noiseless
intermediate images, use ``sim.optical_train.image_planes[0].data``.

---

Instrument packages and YAML files
------------------------------------

All instrument-specific data lives in the IRDB, downloaded via::

    scopesim.download_packages(["Armazones", "ELT", "MORFEO", "MICADO"])

Each package is a directory of YAML files and associated data (FITS tables,
transmission curves, PSF cubes). A ``default.yaml`` defines which YAML files
to load for each observing mode. The YAML files list the ``effects`` that make
up the optical train, with their parameters.

You can browse the downloaded packages in ``./inst_pkgs/`` to understand what
parameters are available and what data files are used.

---

Putting it all together
------------------------

A complete simulation::

    import scopesim
    import scopesim_templates as st
    from scopesim import Simulation

    # Download instrument packages (once)
    scopesim.download_packages(["Armazones", "ELT", "MORFEO", "MICADO"])

    # 1. Create the on-sky target
    src = st.stellar.cluster(mass=1e4, distance=50000, filter_name="V")

    # 2. Load the optical train for the chosen mode
    sim = Simulation("MICADO", mode=["MCAO_4mas", "IMG"])

    # 3. Override any parameters you want to change
    sim.settings["!OBS.dit"] = 120
    sim.settings["!OBS.ndit"] = 5
    sim.settings["!OBS.filter_name"] = "Ks"

    # 4. Run the simulation — returns an astropy HDUList
    hdu = sim(src)
    hdu.writeto("my_simulation.fits", overwrite=True)

    # Optional: inspect the noiseless focal-plane image
    noiseless = sim.optical_train.image_planes[0].data
