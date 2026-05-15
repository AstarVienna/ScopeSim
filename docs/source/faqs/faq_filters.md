# Filters

## How do I list available filters?

```python
sim.optical_train["filter_wheel"].filters
```

This returns a table of filter names and their wavelength coverage. The exact
name of the filter wheel effect may differ per instrument — use
`sim.optical_train.effects` to see what effects are loaded and find the
relevant one.

---

## How do I change the filter?

```python
sim.settings["!OBS.filter_name"] = "Ks"
```

Available filter names are shown by `sim.optical_train["filter_wheel"].filters`.
The parameter path (`!OBS.filter_name`) may differ slightly between instruments
— check `sim.settings["!OBS"]` if the default path does not work.

---

## How do I plot a filter transmission curve?

```python
filt = sim.optical_train["filter_wheel"].current_filter
filt.plot()
```

Or access the transmission table directly:

```python
filt.table   # astropy Table with "wavelength" and "transmission" columns
```

---

## How do I use a custom filter?

Create a `FilterCurve` effect from a two-column ASCII file (wavelength [µm],
transmission [0–1]):

```python
from scopesim.effects import FilterCurve

my_filter = FilterCurve(
    filename="my_filter.dat",
    name="my_custom_filter",
)
sim.optical_train.optics_manager.add_effect(my_filter)
```

Or pass a filter directly as a `synphot.SpectralElement`:

```python
from synphot import SpectralElement, Empirical1D
import astropy.units as u
import numpy as np

wave = np.linspace(2.0, 2.5, 100) * u.um
trans = np.exp(-0.5 * ((wave.value - 2.2) / 0.1) ** 2)
sp = SpectralElement(Empirical1D, points=wave, lookup_table=trans)
```

---

## What filters are available for MICADO?

MICADO's standard filter set includes:

| Name | Wavelength range | Notes |
|---|---|---|
| `J` | 1.17–1.33 µm | Broadband J |
| `H` | 1.49–1.78 µm | Broadband H |
| `Ks` | 2.00–2.37 µm | Broadband Ks |
| `Y` | 0.97–1.07 µm | Y-band |
| `z` | 0.85–0.95 µm | z-band |
| `Br-gamma` | 2.16 µm | Narrow-band Brγ emission |
| `H2` | 2.12 µm | Narrow-band H₂ emission |
| `FeII` | 1.64 µm | Narrow-band [Fe II] emission |

Use `sim.optical_train["filter_wheel"].filters` after loading the MICADO
package for the authoritative current list.

---

## Why does my flux change when I change wavelength range?

ScopeSim integrates over wavelengths. If `!SIM.spectral.wave_min` or
`wave_max` is set narrower than the filter bandpass, flux outside the
simulation range is silently dropped. Always ensure the simulation wavelength
range covers at least the full filter bandpass:

```python
sim.settings["!SIM.spectral.wave_min"] = 1.9   # µm
sim.settings["!SIM.spectral.wave_max"] = 2.5   # µm
```
