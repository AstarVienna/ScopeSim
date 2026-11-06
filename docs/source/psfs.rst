Point Spread Functions
=======================

ScopeSim models the PSF as one or more **PSF effects** in the optical train.
Several PSF classes are available, each suited to different simulation
scenarios. Instrument packages in the IRDB configure the appropriate PSF
automatically, but this page explains the options so you can understand what
is active and how to modify it.

---

Available PSF effect types
---------------------------

FieldConstantPSF
~~~~~~~~~~~~~~~~~

A PSF that is the same across the entire field of view, read from a
FITS file. The FITS file may contain a single kernel or a wavelength-
dependent PSF cube (one kernel per wavelength slice).

**Best for:** Instruments with a relatively uniform PSF, or when you want a
single representative PSF for speed.

**YAML example:**

.. code-block:: yaml

    - name: micado_psf
      class: FieldConstantPSF
      kwargs:
        filename: "PSF_MICADO_SCAO_1.5mas.fits"

FieldVaryingPSF
~~~~~~~~~~~~~~~~

A PSF that varies with position across the field, read from a FITS file
that contains PSF kernels on a spatial grid (and optionally also as a
function of wavelength). ScopeSim interpolates between grid points to
evaluate the PSF at any detector position.

**Best for:** AO PSFs that degrade away from the guide star, such as SCAO
systems where the PSF varies strongly across the MICADO field.

.. note::

   Field-varying PSF files for MICADO (SCAO modes) are large (~1–5 GB) and
   must be downloaded separately from the regular instrument packages::

       scopesim.download_psfs(psf_name="MICADO_SCAO_4mas_FV")

AnisocadoConstPSF
~~~~~~~~~~~~~~~~~~

A field-constant AO PSF generated on-the-fly from the
`AnisoCADO <https://anisocado.readthedocs.io/en/latest/>`_ package.
AnisoCADO computes SCAO PSFs for the ELT using a semi-analytical model
parameterised by atmospheric conditions and guide star separation.

**Best for:** Generating realistic MICADO SCAO PSFs without downloading
large pre-computed files. Useful for exploring how the PSF depends on
observing conditions.

**Requirements:** ``pip install anisocado``

**Example:**

.. code-block:: python

    from scopesim.effects import AnisocadoConstPSF

    psf = AnisocadoConstPSF(
        filename=None,
        strehl_ratio=0.5,
        wavelength=2.15,    # [µm]
        pixel_scale=0.004,  # [arcsec/pix]
    )

SeeingPSF
~~~~~~~~~~

An analytical seeing-limited PSF modelled as a 2D Gaussian with FWHM
read from ``!ATMO.seeing``. Simple and fast — no PSF file needed.

**Best for:** Seeing-limited simulations, or as a quick substitute while
developing a simulation before loading the real AO PSF.

**Example (override the seeing FWHM):**

.. code-block:: python

    sim.settings["!ATMO.seeing"] = 0.8   # [arcsec]

GaussianDiffractionPSF
~~~~~~~~~~~~~~~~~~~~~~~

An analytical PSF combining a Gaussian seeing component with a diffraction-
limited Airy disk. Parameterised by the telescope diameter and atmospheric
seeing.

**Best for:** Quick exploratory simulations or instruments without a dedicated
PSF file.

---

Disabling or replacing the PSF
--------------------------------

To check what PSF effect is loaded::

    sim.optical_train.effects   # look for effects of type *PSF

To disable the PSF (simulate a delta-function PSF)::

    sim.optical_train["psf_effect_name"].include = False

Replace with a simpler PSF for exploratory runs::

    from scopesim.effects import SeeingPSF
    sim.optical_train.optics_manager.add_effect(
        SeeingPSF(fwhm=0.8, name="quick_psf")
    )
    sim.optical_train["original_psf"].include = False

---

Field-varying PSFs for MICADO
-------------------------------

MICADO SCAO observations use a field-varying PSF — the PSF is excellent
near the guide star and degrades toward the field edge. Pre-computed
field-varying PSF cubes (generated with AnisoCADO from ESO SCAO simulations)
are available for download separately from the IRDB package.

.. list-table::
   :widths: 30 30 40
   :header-rows: 1

   * - PSF file
     - Pixel scale
     - Description
   * - ``MICADO_SCAO_4mas_FV``
     - 4 mas/pix
     - 7×7 spatial grid, H+K bands
   * - ``MICADO_SCAO_1.5mas_FV``
     - 1.5 mas/pix
     - 7×7 spatial grid, H+K bands

These files are several GB each. Use ``FieldConstantPSF`` (the package
default) for most simulations; switch to ``FieldVaryingPSF`` when the
spatial variation of the PSF matters for your science case (e.g. astrometric
calibration, PSF-subtracted imaging).

---

PSF file format
----------------

ScopeSim reads PSF FITS files in two formats:

**Single kernel** — a 2D FITS image with:

- ``PIXSCALE`` keyword in the header giving the kernel pixel scale [arcsec]
- Pixel values representing the PSF intensity (will be normalised to sum=1)

**Wavelength cube** — a 3D FITS cube with:

- First axis: wavelength
- Second and third axes: spatial kernel
- ``WAVE0``, ``DWAVE``, ``WAVEUNIT`` header keywords defining the wavelength axis

For field-varying PSFs, the FITS file contains a grid of kernels arranged by
spatial position. See the IRDB documentation for the specific format used by
MICADO field-varying PSF files.
