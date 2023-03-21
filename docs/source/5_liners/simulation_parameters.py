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

# # Global rc simulation parameters
#
# Default global simulation parameters used as a base layer for all instruments. 
#
# Also accessible via an empty `scopesim.UserCommands()` or via an empty `scopesim.OpticalTrain()`object

# +
import scopesim

scopesim.rc.__currsys__["!SIM.random.seed"] = 9001
scopesim.rc.__currsys__["!SIM.file.local_packages_path"] = "./"

scopesim.rc.__currsys__["!SIM"]
