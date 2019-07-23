# ScopeSim 
#### A telescope observation simulator for Python

[![Build Status](https://travis-ci.org/astronomyk/ScopeSim.svg?branch=master)](https://travis-ci.org/astronomyk/ScopeSim)
[![Documentation Status](https://readthedocs.org/projects/scopesim/badge/?version=latest)](https://scopesim.readthedocs.io/en/latest/?badge=latest)

[![Python 2.7](https://img.shields.io/badge/Python-2.7-red.svg)]()
[![Python 3.5](https://img.shields.io/badge/Python-3.5-brightgreen.svg)]()
[![Python 3.6](https://img.shields.io/badge/Python-3.6-brightgreen.svg)]()
[![Python 3.7](https://img.shields.io/badge/Python-3.7-brightgreen.svg)]()


## Summary

Telescopy aims to simulate images of astronomical objects observed with visual 
and infrared instruments. It does this by creating models of the optical train 
and astronomical objects and then pushing the object through the optical train. 
The resulting 2D image is then broadcast to a detector chip and read out into a 
FITS file. 

This code was originally based on the [SimCADO](www.univie.ac.at/simcado) package

## Dependencies

```
numpy >= 1.13
scipy
astropy
synphot
pyyaml
requests
beautifulsoup
```

## Documentation
https://scopesim.readthedocs.io/en/latest/
