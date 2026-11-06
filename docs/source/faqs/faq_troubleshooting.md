# Troubleshooting

Common errors and how to fix them.

---

## `File cannot be found: default.yaml`

**Cause:** ScopeSim cannot find the instrument packages. They have not been
downloaded, or you are running from the wrong directory.

**Fix:** Download the packages into your current working directory:

```python
import scopesim
scopesim.download_packages(["Armazones", "ELT", "MORFEO", "MICADO"])
```

Or tell ScopeSim where the packages live:

```python
scopesim.rc.__config__["!SIM.file.local_packages_path"] = "/path/to/inst_pkgs"
```

---

## `RuntimeError: No package named X found`

**Cause:** The package `X` has not been downloaded or the name is wrong.

**Fix:** Check available packages:

```python
scopesim.list_packages()           # packages available on the server
scopesim.list_packages(local=True) # packages already installed locally
```

---

## Simulation returns all zeros / black image

**Cause 1:** The source is outside the field of view.

**Fix:** Check your source coordinates. The field centre is typically at
`(0, 0)` arcsec unless you have set a WCS offset. Point sources need
`x` and `y` in arcsec relative to the field centre.

**Cause 2:** The source spectrum has no flux in the simulation wavelength
range.

**Fix:** Check `sim.settings["!SIM.spectral.wave_min"]` and `wave_max` match
the filter you are using.

---

## `KeyError: !OBS.some_parameter`

**Cause:** The bang-string path does not exist in the loaded YAML package.

**Fix:** Inspect the settings to find the correct key:

```python
sim.settings["!OBS"]    # view all OBS-level parameters
sim.settings            # view the full settings dict
```

Bang-strings are case-sensitive. Use `!OBS.dit` not `!OBS.DIT`.

---

## Simulation is very slow

**Cause 1:** The spatial oversampling is high. The default chunk size may
create many small FOV tiles.

**Fix:** Increase the chunk size:

```python
sim.settings["!SIM.computing.chunk_size"] = 4096   # default 2048
```

**Cause 2:** A `FieldVaryingPSF` with many spatial positions is loaded.

**Fix:** Substitute a simpler PSF for exploratory runs:

```python
sim.optical_train["psf"].include = False
```

**Cause 3:** Sub-pixel mode is on.

**Fix:**

```python
sim.settings["!SIM.sub_pixel.flag"] = False
```

---

## `AttributeError: 'Simulation' object has no attribute 'X'`

**Cause:** You are calling a SimCADO-style attribute on a ScopeSim object.

**Fix:** See the [SimCADO migration guide](../simcado_migration.rst) for the
ScopeSim equivalent.

---

## Output FITS has unexpected number of extensions

Each chip in the detector array produces one FITS extension. A 3×3 H2RG
array (like MICADO) produces 9 science extensions plus a primary HDU, giving
10 extensions total. Access them by index:

```python
hdu = sim(src)
len(hdu)          # total number of extensions
hdu[1].data       # chip 1 pixel data
hdu[1].header     # chip 1 header with WCS
```

---

## Getting a bug report

When reporting issues, always include:

```python
import scopesim
print(scopesim.bug_report())
```

File issues at https://github.com/AstarVienna/ScopeSim/issues
