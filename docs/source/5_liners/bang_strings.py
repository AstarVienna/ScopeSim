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

# # Using !-string and #-string commands
#
# ## !-strings are for setting simulation parameters
#
# ### TL;DR
#
#     import scopesim as sim
#     opt = sim.load_example_optical_train()
#     opt.cmds["!ATMO"]
#     opt.cmds["!ATMO.background"]
#     opt.cmds["!ATMO.background.filter_name"]
#
# .. note: !-strings only work on `UserCommands` objects
#
# !-strings are a convenient way of accessing multiple layers of a nested dictionary structure with a single string using the format:
#
#     "!<ALIAS>.<sub-dict>...<sub-dict>.<param>"
#     
# Any level of the nested dictionary can be reached by truncating the keyword.

import scopesim as sim
opt = sim.load_example_optical_train()

opt.cmds["!ATMO"]

opt.cmds["!ATMO.background"]

opt.cmds["!ATMO.background.filter_name"]

# ## #-strings are for accessing Effect object parameters
#
# ### TL;DR
#
#     opt.effects
#     opt["#exposure_action."]
#     opt["#exposure_action.ndit"]
#     opt["#exposure_action.ndit!"]
#
#
# .. note: !-strings only work on `OpticalTrain` objects
#
# Similar to !-strings, #-strings allow us to get at the preset values inside the Effect-objects of the optical system. #-strings allow us to pring the contents of an effect's meta dictionary.
#
# First let's list the effects

opt.effects

# We list the meta dictionary contents by using the string format 
#
#     "#<effect-name>."
#     
# .. note: The `.` at the end is important, otherwise the optical train will look for a non-existant effect named `#<effect-name>`

opt["#exposure_action."]

# We print a specific meta parameter by adding it after the `.`

opt["#exposure_action.ndit"]

# Notice that the value of this dictionary entry is itself a !-string. We can resolve this by adding a `!` to the end of the string, to force it to get the actual value from `opt.cmds`:

opt["#exposure_action.ndit!"]
