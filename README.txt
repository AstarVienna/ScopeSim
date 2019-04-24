# ScopeSim 
#### A telescope observation simulator for Python

[![Build Status](https://travis-ci.org/astronomyk/TelescoPy.svg?branch=master)](https://travis-ci.org/astronomyk/TelescoPy.svg?branch=master)
[![Documentation Status](https://readthedocs.org/projects/telescopy/badge/?version=latest)](https://telescopy.readthedocs.io/en/latest/?badge=latest)

## Summary

Telescopy aims to simulate images of astronomical objects observed with visual 
and infrared instruments. It does this by creating models of the optical train 
and astronomical objects and then pushing the object through the optical train. 
The resulting 2D image is then broadcast to a detector chip and read out into a 
FITS file. 

This code was originally based on the [SimCADO](www.univie.ac.at/simcado) package

## Dependencies

```
numpy
scipy
astropy
synphot
requests
```

## Documentation
https://telescopy.readthedocs.io/en/latest/
