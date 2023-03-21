# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Turning Effect objects on or off
#
# **TL;DR**
#
#     optical_train = sim.load_example_optical_train()
#     
#     optical_train.effects
#     optical_train["detector_linearity"].include = False
#     optical_train["detector_linearity"].meta["include"] = True
#
#
# To list all the effects in an optical train, we do use the `effects` attribute.
#
# Alternatively, we can call `opt.optics_manager.all_effects()`

# +
import scopesim as sim

opt = sim.load_example_optical_train()
opt.effects
# -

# Turning an effect on or off is as simple as using setting the `.include` attribute to `true` or `False`:

opt["slit_wheel"].include = True
opt["slit_wheel"].include = False
