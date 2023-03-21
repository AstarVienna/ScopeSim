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

# # Getting started
#
# <div class="alert alert-info">
#
# Note: For instrument specific guides, please see the [IRDB](https://irdb.readthedocs.io/en/latest/)
#
# </div>
#
# A basic simulation would look something like this:

# +
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import scopesim as sim
from scopesim.source import source_templates as st

src = st.star_field(n=100, 
                    mmax=15,      # [mag]
                    mmin=20, 
                    width=200)    # [arcsec]

opt = sim.load_example_optical_train()
opt.cmds["!OBS.dit"] = 60         # [s]
opt.cmds["!OBS.ndit"] = 10

opt.observe(src)
hdulist = opt.readout()[0]

plt.figure(figsize=(10,8))
plt.imshow(hdulist[1].data, norm=LogNorm(vmin=1))
plt.colorbar()
# -

# ## Code breakdown
#
# Let's break this down a bit.
#
# There are three major components of any simulation workflow:
#
# 1. the target description,
# 2. the telescope/instrument model, and
# 3. the observation.
#
# For the target description we are using the ScopeSim internal template functions from `scopesim.source.source_templates`, however many more dedicated science related templates are available in the external python package [ScopeSim-Templates](https://scopesim-templates.readthedocs.io/en/latest/)
#
# Here we create a field of 100 A0V stars with Vega magnitudes between V=15 and V=20 within a box of 200 arcsec:

src = st.star_field(n=100, 
                    mmax=15,      # [mag]
                    mmin=20, 
                    width=200)    # [arcsec]

# Next we load the sample optical train object from ScopeSim.
#
# Normally we will want to use an actual instrument. Dedicated documentation for real telescope+instrument systems can be found in the documentation sections of the individual instruments in the [Instrument Reference Database (IRDB) documentation](https://irdb.readthedocs.io/en/latest/)
#
# For real instruments loading the optical system generally follows a different pattern:
#
#     cmd = sim.UserCommands(use_instrument="instrument_name", set_modes=["mode_1", "mode_2"])
#     opt = sim.OpticalTrain(cmds)
#
# Once we have loaded the instrument, we can set the observation parameters by accessing the internal commands dictionary:

opt = sim.load_example_optical_train(set_modes=["imaging"])
opt.cmds["!OBS.dit"] = 60         # [s]
opt.cmds["!OBS.ndit"] = 10

# Finally we observe the target source and readout the detectors.
#
# What is returned (`hdulist`) is an `astropy.fits.HDUList` object which can be saved to disk in the standard way, or manipulated in a python session.

opt.observe(src)
hdulist = opt.readout()[0]

# ## Tips and tricks
#
# ### Focal plane images
#
# Intermediate frames of the focal plane image without the noise proerties can be accessed by looking inside the optical train object and accessing the first image plane:

noiseless_image = opt.image_planes[0].data

# ### Turning optical effects on or off
#
# All effects modelled by the optical train can be listed with the `.effects` attribute:

opt.effects

# These can be turned on or off by using their name and the `.include` attribute:

opt["detector_linearity"].include = False

# ### Listing available modes and filters
#
# The list of observing modes can be found by using the `.modes` attribute of the commands objects:

opt.cmds.modes

# The names of included filters can be found in the filter effect. Use the name of the filter object from the table above to list these:

opt["filter_wheel"].filters

# ### Setting observation sequence
#
# Although this could be different for some instruments, most instruments use the `exptime = ndit * dit` format.
# `ndit`and `dit` are generally accessible in the top level `!OBS` dictionary of the command object in the optical train.

opt.cmds["!OBS.dit"] = 60         # [s]
opt.cmds["!OBS.ndit"] = 10

# ### Listing and changing simulation parameters
#
# The command dictionary inside the optical system contains all the necessary paramters.

opt.cmds

# The command object is a series of nested dictionaries that can be accessed using the `!-string` format:
#
#     opt.cmds["!<alias>.<param>"]
#     opt.cmds["!<alias>.<sub_dict>.<param>"]
#     
# For example, setting the atmospheric background level is achieved thusly:

opt.cmds["!ATMO.background.filter_name"] = "K"
opt.cmds["!ATMO.background.value"] = 13.6

# ## More information
#
# For more information on how to use ScopeSim be see:
#
# - [Use Examples](examples/index.rst)
# - [Instrument Specific Documentation](https://irdb.readthedocs.io/en/latest/)
# - [Effect data formats](effects/formats)
# - [Setting up a custom instrument]()
#
#
# ## Contact
#
# - For bugs, please add an [issue to the github repo](https://github.com/AstarVienna/ScopeSim/issues)
# - For enquiries on implementing your own instrument package, please drop us a line at astar.astro@univie.ac.at or kieran.leschinski@univie.ac.at
