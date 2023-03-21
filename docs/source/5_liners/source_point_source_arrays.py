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

# # Source : Point sources from arrays
#
# Collections of point sources can be initialised through either a collection of arrays, or an astropy Table

# +
import numpy as np
import astropy.table as table
from astropy import units as u

import scopesim

# how many stars
n = 200

# random coordinated in a 100 arcsec box
x, y  = 100 * np.random.random(size=(2, n)) - 50

# All stars reference the Vega spectrum
ref = np.zeros(n)      
# Each star is given a different weighting, i.e. magnitude
weight = 10**(-0.4*np.linspace(10, 20, n))

# Note: The Pyckles and SpeXtra libraries contain many more stellar and galactic spectra
vega = scopesim.source.source_templates.vega_spectrum(mag=20)
# -

# ## astropy.Table

# +
tbl = table.Table(names=["x", "y", "ref", "weight"],
                  data= [x,    y,   ref,   weight],
                  units=[u.arcsec, u.arcsec, None, None])

table_source = scopesim.Source(table=tbl, spectra=[vega])

table_source.plot()
# -

# ## From loose arrays

# +
point_source = scopesim.Source(spectra=[vega], x=x, y=y, ref=ref, weight=weight)

point_source.plot()
