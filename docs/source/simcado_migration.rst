Migrating from SimCADO
======================

ScopeSim is the successor to SimCADO. While the underlying simulation
philosophy is the same, the architecture has been redesigned to be instrument-
agnostic. This page maps the old SimCADO API onto ScopeSim equivalents.

.. note::

   SimCADO only simulated MICADO. ScopeSim is a general framework — to
   simulate MICADO you need the ScopeSim engine plus the MICADO instrument
   package from the IRDB.

---

The minimum working example
----------------------------

**SimCADO** ::

    import simcado
    src = simcado.source.star_field(100, 15, 20, width=10)
    simcado.run(src, filename="output.fits")

**ScopeSim** ::

    import scopesim
    import scopesim_templates as st
    from scopesim import Simulation

    scopesim.download_packages(["Armazones", "ELT", "MORFEO", "MICADO"])

    src = st.stellar.star_field(100, "V", 15, 20, width=10)
    sim = Simulation("MICADO", mode=["MCAO_4mas", "IMG"])
    sim(src, dit=60, ndit=10)


---

API mapping
-----------

Installation / setup
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - SimCADO
     - ScopeSim equivalent
   * - ``pip install simcado``
     - ``pip install scopesim scopesim_templates``
   * - ``simcado.get_extras()``
     - ``scopesim.download_packages(["Armazones", "ELT", "MORFEO", "MICADO"])``
   * - Data files bundled with package
     - Instrument packages live in ``./inst_pkgs/`` after download

Source creation
~~~~~~~~~~~~~~~

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - SimCADO
     - ScopeSim equivalent
   * - ``simcado.source.star(mag, filter_name)``
     - ``scopesim_templates.stellar.star(mag, filter_name, ...)``
   * - ``simcado.source.stars(mags, x, y, filter_name)``
     - ``scopesim_templates.stellar.stars(filter_name, mags, ...)``
   * - ``simcado.source.star_grid(n, mag_min, mag_max)``
     - ``scopesim_templates.stellar.star_grid(n, mag_min, mag_max, ...)``
   * - ``simcado.source.cluster(mass, distance)``
     - ``scopesim_templates.stellar.cluster(mass, distance, ...)``
   * - ``simcado.source.source_from_image(imgs, lam, spectra)``
     - ``scopesim.Source(image_hdu=fits_image, spectra=spectra, ...)``
   * - ``src1 + src2``
     - ``src1 + src2`` (unchanged)

Configuration / commands
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - SimCADO
     - ScopeSim equivalent
   * - ``cmds = simcado.UserCommands()``
     - ``sim.settings`` (accessed from the ``Simulation`` object)
   * - ``cmds["OBS_EXPTIME"] = 3600``
     - ``sim.settings["!OBS.dit"] = 360; sim.settings["!OBS.ndit"] = 10``
   * - ``cmds["OBS_DIT"] = 60``
     - ``sim.settings["!OBS.dit"] = 60``
   * - ``cmds["OBS_NDIT"] = 10``
     - ``sim.settings["!OBS.ndit"] = 10``
   * - ``cmds["INST_FILTER_TC"] = "K"``
     - ``sim.settings["!OBS.filter_name"] = "K"``
   * - ``cmds["SIM_PIXEL_SCALE"] = 0.004``
     - Determined by the chosen observing mode (e.g. ``"MCAO_4mas"``)
   * - ``cmds["ATMO_SEEING"] = 0.8``
     - ``sim.settings["!ATMO.seeing"] = 0.8``
   * - Keyword names like ``OBS_EXPTIME``
     - Bang-string notation: ``!OBS.dit``, ``!SIM.spectral.wave_min``, etc.

Running simulations
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - SimCADO
     - ScopeSim equivalent
   * - ``simcado.run(src, cmds=cmds, filename="out.fits")``
     - ``sim(src)`` — returns an ``astropy.fits.HDUList``
   * - ``ot = simcado.OpticalTrain(cmds); ot.observe(src)``
     - ``sim.optical_train.observe(src); sim.optical_train.readout()``
   * - ``det = simcado.Detector(cmds); det.read_out()``
     - Handled internally by ``sim(src)``
   * - ``det.array`` (raw pixel data)
     - ``sim(src)[1].data`` (FITS extension 1)

Observing modes
~~~~~~~~~~~~~~~

SimCADO was MICADO-specific and configured via flat keywords. ScopeSim uses
named observing modes defined in the YAML package files.

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - SimCADO configuration
     - ScopeSim equivalent
   * - ``cmds["SIM_PIXEL_SCALE"] = 0.004`` (4 mas)
     - ``Simulation("MICADO", mode=["MCAO_4mas", "IMG"])``
   * - ``cmds["SIM_PIXEL_SCALE"] = 0.0015`` (1.5 mas)
     - ``Simulation("MICADO", mode=["SCAO_1.5mas", "IMG"])``
   * - ``cmds["SCOPE_USE_MIRROR_BG"] = True``
     - Controlled by an ``SurfaceList`` effect in the YAML package

Listing available options::

    sim.settings.modes          # available observing modes
    sim.optical_train["filter_wheel"].filters   # available filters
    sim.optical_train.effects                   # all active effects

Optical train inspection
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - SimCADO
     - ScopeSim equivalent
   * - ``ot.apply_optical_train(src, det)``
     - ``sim.optical_train.observe(src)``
   * - Accessing intermediate images: not straightforward
     - ``sim.optical_train.image_planes[0].data``
   * - Turn effect off: not directly supported
     - ``sim.optical_train["effect_name"].include = False``
   * - List effects: not supported
     - ``sim.optical_train.effects``

---

Key conceptual changes
-----------------------

Modular instrument packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SimCADO was self-contained — all MICADO data was shipped with the Python
package. ScopeSim separates the simulation engine from the instrument data:

- **ScopeSim** — the simulation engine (this package)
- **IRDB** — instrument packages downloaded on first use
- **ScopeSim Templates** — helper functions for creating on-sky sources

This means you must call ``scopesim.download_packages(...)`` before running any
real instrument simulation. The packages are cached in ``./inst_pkgs/`` and
only need to be downloaded once.

Bang-string parameter access
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SimCADO used flat uppercase keyword names (``OBS_EXPTIME``, ``INST_FILTER_TC``).
ScopeSim organises parameters in a nested YAML hierarchy accessed with
bang-string notation::

    sim.settings["!OBS.dit"]              # observation DIT
    sim.settings["!OBS.filter_name"]      # filter selection
    sim.settings["!ATMO.seeing"]          # atmospheric seeing
    sim.settings["!SIM.spectral.wave_min"] # simulation wavelength range

The ``!`` prefix signals a YAML-path lookup. Keys are structured as
``!<ALIAS>.<section>.<parameter>``.

MICADO now needs MORFEO
~~~~~~~~~~~~~~~~~~~~~~~~

In SimCADO, MICADO's adaptive optics (MAORY at the time) was implicit. In
ScopeSim, MICADO works with MORFEO as a separate instrument package that must
be downloaded explicitly::

    scopesim.download_packages(["Armazones", "ELT", "MORFEO", "MICADO"])

For purely diffraction-limited or seeing-limited simulations (no AO), use the
standalone MICADO modes.
