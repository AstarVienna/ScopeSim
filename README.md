# ScopeSim 
## A telescope observation simulator for Python

[![Build Status](https://github.com/AstarVienna/ScopeSim/actions/workflows/tests.yml/badge.svg)](https://github.com/AstarVienna/ScopeSim/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/scopesim/badge/?version=latest)](https://scopesim.readthedocs.io/en/latest/?badge=latest)

[![Build Status](http://github-actions.40ants.com/AstarVienna/ScopeSim/matrix.svg)](https://github.com/AstarVienna/ScopeSim)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Summary

Telescopy aims to simulate images of astronomical objects observed with visual 
and infrared instruments. It does this by creating models of the optical train 
and astronomical objects and then pushing the object through the optical train. 
The resulting 2D image is then broadcast to a detector chip and read out into a 
FITS file. 

This code was originally based on the [SimCADO](https://github.com/astronomyk/simcado) package

## Documentation
The main set of documentation can be found here: 
https://scopesim.readthedocs.io/en/latest/

A basic Jupyter Notebook can be found here: 
[scopesim_basic_intro.ipynb](docs/source/examples/1_scopesim_intro.ipynb)
