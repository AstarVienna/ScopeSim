# ScopeSim 
## A telescope observation simulator for Python

[![Build Status](https://github.com/AstarVienna/ScopeSim/actions/workflows/tests.yml/badge.svg)](https://github.com/AstarVienna/ScopeSim/actions/workflows/tests.yml/badge.svg)
[![Build Status](https://github.com/AstarVienna/ScopeSim/actions/workflows/notebooks_with_irdb_download.yml/badge.svg)](https://github.com/AstarVienna/ScopeSim/actions/workflows/notebooks_with_irdb_download.yml/badge.svg)
[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
![dev version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FAstarVienna%2FScopeSim%2Fmain%2Fpyproject.toml&query=%24.tool.poetry.version&label=dev%20version&color=teal)

[![Documentation Status](https://readthedocs.org/projects/scopesim/badge/?version=latest)](https://scopesim.readthedocs.io/en/latest)
[![codecov](https://codecov.io/gh/AstarVienna/ScopeSim/graph/badge.svg)](https://codecov.io/gh/AstarVienna/ScopeSim)
[![PyPI - Version](https://img.shields.io/pypi/v/ScopeSim)](https://pypi.org/project/ScopeSim/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ScopeSim)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Citation](https://img.shields.io/badge/DOI-10.1117%2F12.2559784-blue)](https://doi.org/10.1117/12.2559784)

## Summary

ScopeSim aims to simulate images of astronomical objects observed with visual 
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
