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

# # 1: A quick use case for MICADO at the ELT
#
#
# ## A brief introduction into using ScopeSim to observe a cluster in the LMC

# +
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# %matplotlib inline

import scopesim as sim
import scopesim_templates as sim_tp

# [Required for Readthedocs] Comment out this line if running locally
tmpdir = TemporaryDirectory()
sim.rc.__config__["!SIM.file.local_packages_path"] = tmpdir.name
# -

# Download the required instrument packages for an observation with MICADO at the ELT

sim.download_packages(["Armazones", "ELT", "MAORY", "MICADO"])

# Create a star cluster using the ``scopesim_templates`` package

cluster = sim_tp.stellar.clusters.cluster(mass=1000,         # Msun
                                          distance=50000,    # parsec
                                          core_radius=0.3,     # parsec
                                          seed=9002)

# Make the MICADO optical system model with ``OpticalTrain``. Observe the cluster ``Source`` object with the ``.observe()`` method and read out the MICADO detectors with ``.readout()``. 
#
# The resulting FITS file can either be returned as an ``astropy.fits.HDUList`` object, or saved to disk using the optional ``filename`` parameter

micado = sim.OpticalTrain("MICADO")
micado.observe(cluster)
hdus = micado.readout()
# micado.readout(filename="TEST.fits")

# Display the contents the first HDU

plt.figure(figsize=(10,8))
plt.imshow(hdus[0][1].data, norm=LogNorm(vmax=3E4, vmin=3E3), cmap="hot")
plt.colorbar()

# ## TL;DR
#
# ```
# import scopesim as sim
# import scopesim_templates as sim_tp
#
# sim.download_packages(["Armazones", "ELT", "MAORY", "MICADO"])
#
# cluster = sim_tp.stellar.clusters.cluster(mass=1000,         # Msun
#                                           distance=50000,    # parsec
#                                           core_radius=0.3,     # parsec
#                                           seed=9002)
#
# micado = sim.OpticalTrain("MICADO")
# micado.observe(cluster)
#
# hdus = micado.readout()
# # micado.readout(filename="TEST.fits")
# ```
