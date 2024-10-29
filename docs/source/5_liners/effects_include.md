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

# Turning Effect objects on or off

To list all the effects in an optical train, we do use the `effects` attribute.

Alternatively, we can call `opt.optics_manager.all_effects()`

```{code-cell} ipython3
import scopesim as sim

opt = sim.load_example_optical_train()
opt.effects
```

Turning an effect on or off is as simple as using setting the `.include` attribute to `true` or `False`:

```{code-cell} ipython3
opt["slit_wheel"].include = True
opt["slit_wheel"].include = False
```
