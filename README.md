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

This code was originally based on the [SimCADO](www.univie.ac.at/simcado) package

## Documentation
The main set of documentation can be found here: 
https://scopesim.readthedocs.io/en/latest/

A basic Jupyter Notebook can be found here: 
[scopesim_basic_intro.ipynb](docs/source/_static/scopesim_basic_intro.ipynb)


## Dependencies

For [![Python 3.6](https://img.shields.io/badge/Python-3.6-brightgreen.svg)]() and above the latest versions of these packages are compatible with ScopeSim:

    numpy >= 1.16
    scipy >= 1.0.0
    astropy >= 2.0
    pyyaml >= 5.1
    requests >= 2.20
    beautifulsoup4 >= 4.4
    synphot >= 0.1.3

For [![Python 3.5](https://img.shields.io/badge/Python-3.5-yellow.svg)]() the following packages may not exceed these version numbers:

    astropy <= 3.2.3
    synphot <= 0.1.3

#### Oldest currently tested system 

[![Python 3.5](https://img.shields.io/badge/Python-3.5-yellow.svg)]()

[![Numpy](https://img.shields.io/badge/Numpy-1.16-brightgreen.svg)]()
[![Astropy](https://img.shields.io/badge/Astropy-2.0-brightgreen.svg)]()
[![Scipy](https://img.shields.io/badge/Scipy-1.0.0-brightgreen.svg)]()

[![Synphot](https://img.shields.io/badge/Synphot-0.1.3-brightgreen.svg)]()
[![requests](https://img.shields.io/badge/requests-2.20.0-brightgreen.svg)]()
[![beautifulsoup4](https://img.shields.io/badge/beautifulsoup4-4.4-brightgreen.svg)]()
[![pyyaml](https://img.shields.io/badge/pyyaml-5.1-brightgreen.svg)]()

#### Things to watch out for with Synphot
Numpy>=1.16 must be used for synphot to work
For Astropy<4.0, only Synphot<=0.1.3 works

#### Optional dependencies
[![skycalc_ipy](https://img.shields.io/badge/skycalc_ipy->=0.1-brightgreen.svg)]()
[![anisocado](https://img.shields.io/badge/anisocado->=0.1-brightgreen.svg)]()
[![Matplotlib](https://img.shields.io/badge/Matplotlib->=1.5-brightgreen.svg)]()

