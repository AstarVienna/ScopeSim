---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Effects Overview

In ScopeSim, **Effects** are the building blocks of an optical system simulation.
Each `Effect` object represents a single physical phenomenon — atmospheric
seeing, mirror reflectivity, filter transmission, detector noise, and so on.
An `OpticalTrain` collects all the effects from an instrument package and
applies them in sequence to transform a `Source` (the on-sky light distribution)
into a realistic detector readout.

Effects are typically defined in YAML instrument configuration files, but can
also be created and added programmatically. Each effect implements an
`apply_to(obj)` method that receives an object at a specific stage of the
simulation pipeline and returns the modified object.

## Listing Effects in an Optical Train

To see all effects loaded in an optical train, use the `.effects` attribute:

```{code-cell} ipython3
import scopesim as sim

opt = sim.load_example_optical_train()
opt.effects
```

## The Simulation Pipeline

ScopeSim processes effects in a strict order determined by each effect's
**z_order** — a numerical priority that maps to a specific pipeline stage.
Effects with lower z_order values run first. Each effect's z_order can contain
multiple values, allowing it to participate in more than one stage (e.g., setup
and application).

The pipeline has three main phases:

```{mermaid}
%%{init: {"theme": "dark"} }%%
flowchart TB
    subgraph Setup ["Setup Phase"]
        S1["FOV Setup\nz = 200..299"]
        S2["Image Plane Setup\nz = 300..399"]
        S3["Detector Setup\nz = 400..499"]
    end
    subgraph Observe ["observe() Phase"]
        O1["Source Effects\nz = 500..599\n<i>TER curves, filters</i>"]
        O2["FOV Effects\nz = 600..699\n<i>PSFs, spectral traces, shifts</i>"]
        O3["Image Plane Effects\nz = 700..799\n<i>Vibration, flat fields</i>"]
    end
    subgraph Readout ["readout() Phase"]
        R1["Detector Effects\nz = 800..899\n<i>Noise, dark current, QE</i>"]
        R2["Detector Array Effects\nz = 900..999\n<i>Exposure integration</i>"]
        R3["FITS Header Effects\nz = 1000+"]
    end
    S1 --> S2 --> S3 --> O1 --> O2 --> O3 --> R1 --> R2 --> R3
```

### Z-Order Reference

| Z-Order Range | Pipeline Stage | Applied To | OpticsManager Property |
|:---:|---|---|---|
| 200–299 | FOV setup | `FovVolumeList` | `fov_setup_effects` |
| 300–399 | Image plane setup | `FovVolumeList` | `image_plane_setup_effects` |
| 400–499 | Detector setup | `FovVolumeList` | `detector_setup_effects` |
| 500–599 | Source effects | `Source` | `source_effects` |
| 600–699 | FOV effects | `FieldOfView` | `fov_effects` |
| 700–799 | Image plane effects | `ImagePlane` | `image_plane_effects` |
| 800–899 | Detector effects | `Detector` | `detector_effects` |
| 900–999 | Detector array effects | `Detector` | `detector_array_effects` |
| 1000+ | FITS header effects | `HDUList` | `fits_header_effects` |

## Effect Categories

### Transmission, Emission, and Reflection (TER) Curves

TER curves describe how optical surfaces transmit, emit, and reflect light as
a function of wavelength. They are the most common type of effect and model
everything from mirror coatings to atmospheric transmission to detector quantum
efficiency.

| Class | Description |
|---|---|
| `TERCurve` | Base wrapper for spectral transmission/emission/reflection curves |
| `SurfaceList` | Combines multiple optical surfaces into a system-level TER curve |
| `AtmosphericTERCurve` | Atmospheric transmission and emission |
| `SkycalcTERCurve` | Atmospheric TER from ESO's SkyCalc service |
| `QuantumEfficiencyCurve` | Detector quantum efficiency vs wavelength |
| `FilterCurve` | Bandpass filter transmission from file |
| `TopHatFilterCurve` | Rectangular (top-hat) filter transmission |
| `FilterWheel` | A wheel of selectable filters |
| `TopHatFilterWheel` | A wheel of selectable top-hat filters |
| `SpanishVOFilterCurve` | Filters from the Spanish Virtual Observatory |
| `DownloadableFilterCurve` | Filters downloadable from remote servers |
| `PupilTransmission` | Pupil plane transmission curves |
| `PupilMaskWheel` | Wheel of pupil masks |
| `ADCWheel` | Atmospheric Dispersion Corrector wheel |

### Apertures and Field Masks

Aperture effects define the on-sky field geometry — imaging windows, spectrograph
slits, and IFU fields.

| Class | Description |
|---|---|
| `ApertureMask` | Defines on-sky window coordinates (imaging, slit, IFU, MOS) |
| `RectangularApertureMask` | Rectangular aperture variant |
| `ApertureList` | Container for multiple apertures |
| `SlitWheel` | A wheel of selectable slits |

### Point Spread Functions (PSFs)

PSF effects model the spatial blurring of point sources due to diffraction,
atmospheric seeing, and optical aberrations.

| Class | Description |
|---|---|
| `Vibration` | Wavelength-independent Gaussian vibration kernel |
| `SeeingPSF` | Atmospheric seeing as a Gaussian kernel |
| `GaussianDiffractionPSF` | Diffraction-limited PSF (Gaussian approximation) |
| `NonCommonPathAberration` | NCPA PSF from wavefront error maps |
| `AnisocadoConstPSF` | SCAO PSF from AnisoCADO at a given Strehl ratio |
| `FieldConstantPSF` | Field-constant PSF loaded from a FITS file |
| `FieldVaryingPSF` | Field-varying PSF loaded from a FITS file |

### Shifts and Atmospheric Dispersion

Shift effects apply wavelength-dependent positional offsets to the light
distribution.

| Class | Description |
|---|---|
| `Shift3D` | Base class for wavelength-dependent positional shifts |
| `AtmosphericDispersion` | Wavelength-dependent atmospheric refraction |
| `AtmosphericDispersionCorrection` | ADC correction for atmospheric dispersion |

### Spectral Traces

Spectral trace effects map 3D spectral data cubes onto the 2D detector plane
for spectrographic modes.

| Class | Description |
|---|---|
| `SpectralTraceList` | Maps spectral cubes to detector plane via trace geometries |
| `SpectralTraceListWheel` | A wheel of selectable spectral trace configurations |
| `SpectralEfficiency` | Grating/dispersion efficiency (blaze function) |

### Detector Geometry

Detector geometry effects define the physical layout of detector chips on
the focal plane.

| Class | Description |
|---|---|
| `DetectorList` | Detector chip positions, sizes, pixel scale, and rotation |
| `DetectorWindow` | A sub-region readout window on a detector |
| `DetectorList3D` | 3D detector array definition for spectroscopic modes |

### Electronic and Noise Effects

These effects model the detector electronics and noise sources that affect the
final readout.

| Class | Description |
|---|---|
| `ShotNoise` | Poissonian photon noise |
| `DarkCurrent` | Thermal dark current |
| `Bias` | Constant bias level |
| `BasicReadoutNoise` | Generic readout noise |
| `PoorMansHxRGReadoutNoise` | Simplified HAWAII detector readout noise model |
| `LinearityCurve` | Detector linearity and saturation |
| `ADConversion` | Analog-to-digital conversion (electrons to ADU) |
| `PixelResponseNonUniformity` | Per-pixel gain variations (flat field) |
| `InterPixelCapacitance` | Inter-pixel capacitance crosstalk kernel |

### Exposure and Readout

Exposure effects handle integration time, auto-exposure, and readout formatting.

| Class | Description |
|---|---|
| `AutoExposure` | Auto-determines DIT/NDIT from saturation limits |
| `ExposureIntegration` | Integrates signal over the exposure time |
| `ExposureOutput` | Formats the readout output |
| `DetectorModePropertiesSetter` | Sets mode-specific detector parameters (MINDIT, FULL_WELL, RON) |

### Detector Pixel Effects

| Class | Description |
|---|---|
| `BinnedImage` | Equal pixel binning |
| `UnequalBinnedImage` | Non-uniform pixel binning |
| `ReferencePixelBorder` | Masks reference pixels at detector edges |
| `Rotate90CCD` | Rotates CCD by integer multiples of 90 degrees |

### Observing Strategies

| Class | Description |
|---|---|
| `ChopNodCombiner` | Combines 4 chop/nod positions (AA, AB, BA, BB) |

### FITS Header Effects

| Class | Description |
|---|---|
| `ExtraFitsKeywords` | Adds custom FITS keywords to output headers |
| `EffectsMetaKeywords` | Adds effect metadata to FITS headers |
| `SourceDescriptionFitsKeywords` | Adds source description keywords |
| `SimulationConfigFitsKeywords` | Adds simulation configuration keywords |

### Other

| Class | Description |
|---|---|
| `Shutter` | Simulates a closed shutter (zeros all pixels) |

## YAML Configuration

Effects are typically defined in YAML instrument packages. Each effect entry
specifies the class name and configuration parameters:

```yaml
effects:
  - name: detector_qe_curve
    description: Quantum efficiency of the battery of detectors
    class: QuantumEfficiencyCurve
    kwargs:
      filename: QE_detector_H2RG.dat

  - name: dark_current
    description: Detector dark current
    class: DarkCurrent
    kwargs:
      value: 0.1       # electrons/s/pixel

  - name: filter_wheel
    class: FilterWheel
    kwargs:
      current_filter: "!OBS.filter_name"
      filter_names: [J, H, Ks]
      filename_format: "filters/TC_filter_{}.dat"
```

Parameters prefixed with `!` (called **bang strings**) are resolved dynamically
from the simulation configuration at runtime. For example, `!OBS.filter_name`
reads the current filter selection from the observation commands.

## Interacting with Effects at Runtime

Effects can be accessed, toggled, and modified after the optical train is loaded.

### Enabling and disabling effects

```{code-cell} ipython3
# Turn off an effect
opt["detector_linearity"].include = False
print("detector_linearity included:", opt["detector_linearity"].include)

# Turn it back on
opt["detector_linearity"].include = True
```

### Inspecting effect metadata

```{code-cell} ipython3
opt["dark_current"].meta
```

### Modifying parameters

```{code-cell} ipython3
# Change the dark current value
opt["dark_current"].meta["value"] = 0.5
print("New dark current:", opt["dark_current"].meta["value"])
```

For more tips on interacting with effects, see:
- [Turning Effects on or off](../5_liners/effects_include.md)
- [Using bang strings and hash strings](../5_liners/bang_strings.md)

## See Also

- [Creating Custom Effects](custom_effects.md) — how to write your own Effect subclasses
- [Custom Effects Example Notebook](../examples/3_custom_effects.ipynb) — a worked example using MICADO and a PointSourceJitter effect
- The auto-generated [API Reference for scopesim.effects](../_autosummary/scopesim.effects.html)
