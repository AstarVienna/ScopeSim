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

# # Downloading packages
#
# .. note: Instrument packages are kept in a separate repository: [the Instrument Reference Database (IRDB)]((https://github.com/astronomyk/irdb))
#
# Before simulating anything we need to get the relevant instrument packages. Packages are split into the following categories
#
# - Locations (e.g. Armazones, LaPalma)
# - Telescopes (e.g. ELT, GTC)
# - Instruments (e.g. MICADO, METIS, MAORY, OSIRIS, MAAT)
#
# We need to amke sure we have all the packages required to built the optical system. E.g. observing with MICADO is useless without including the ELT.

# +
import scopesim as sim

from tempfile import TemporaryDirectory
tmpdir = TemporaryDirectory()
sim.rc.__config__["!SIM.file.local_packages_path"] = tmpdir.name
# -

# ## scopesim.download_packages()
#
# ### Stable packages
#
# The simplest way is to simply get the latest stable versions of the packages by calling their names.
#
# Call `list_packages()` or see the [IRDB]((https://github.com/astronomyk/irdb)) for names.

sim.list_packages()

sim.download_packages(["Armazones", "ELT", "MICADO"])

# ### Development version
#
# Use the `release="latest"` parameter to get the latest stable development verions of the packages

sim.download_packages("test_package", release="latest")

# ### Bleeding-edge versions from GitHub
#
# We can also pull the latest verisons of a package from GitHub. This is useful if you are the one writing the package and want to test how it works on a different machine.
#
# The following strings will work:
#
# - Latest from a branch: `release="github:<branch-name>"`
# - From a specific tag: `release="github:<tag>"`
# - From a specific commit: `release="github:<super-long-commit-hash>"`

sim.download_packages("test_package", release="github:dev_master")

sim.download_packages("LFOA", release="github:3c136cd59ceeca551c01c6fa79f87377997f33f9")
