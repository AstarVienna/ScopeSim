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

# Using !-string and #-string commands

## !-strings are for setting simulation parameters

!-strings are a convenient way of accessing multiple layers of a nested dictionary structure with a single string using the format:

    "!<ALIAS>.<sub-dict>...<sub-dict>.<param>"

Any level of the nested dictionary can be reached by truncating the keyword.

**Note: !-strings only work on `UserCommands` objects**

Below is an example of how to use !-strings, using the example optical train.

```{code-cell} ipython3
import scopesim as sim
opt = sim.load_example_optical_train()
```

```{code-cell} ipython3
opt.cmds["!ATMO"]
```

```{code-cell} ipython3
opt.cmds["!ATMO.background"]
```

```{code-cell} ipython3
opt.cmds["!ATMO.background.filter_name"]
```

## #-strings are for accessing Effect object parameters

Similar to !-strings, #-strings allow us to get at the preset values inside the Effect-objects of the optical system. #-strings allow us to pring the contents of an effect's meta dictionary.

**Note: !-strings only work on `OpticalTrain` objects**

Here, we're again using the example optical train defined above. First let's list the effects:

```{code-cell} ipython3
opt.effects
```

We list the meta dictionary contents by using the string format

    "#<effect-name>."

**Note: The `.` at the end is important, otherwise the optical train will look for a non-existant effect named `#<effect-name>`**

```{code-cell} ipython3
opt["#exposure_integration."]
```

We print a specific meta parameter by adding it after the `.`

```{code-cell} ipython3
opt["#exposure_integration.ndit"]
```

Notice that the value of this dictionary entry is itself a !-string. We can resolve this by adding a `!` to the end of the string, to force it to get the actual value from `opt.cmds`:

```{code-cell} ipython3
opt["#exposure_integration.ndit!"]
```
