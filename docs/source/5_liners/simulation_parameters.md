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

# Global rc simulation parameters

Default global simulation parameters used as a base layer for all instruments. 

Also accessible via an empty `scopesim.UserCommands()` or via an empty `scopesim.OpticalTrain()`object

```{code-cell} ipython3
import scopesim

scopesim.rc.__currsys__["!SIM.random.seed"] = 9001
scopesim.rc.__currsys__["!SIM.file.local_packages_path"] = "./"

scopesim.rc.__currsys__["!SIM"]
```
