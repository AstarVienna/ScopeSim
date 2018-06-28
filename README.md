# TelescoPy

[![Build Status](https://travis-ci.org/astronomyk/TelescoPy.svg?branch=master)](https://travis-ci.org/astronomyk/TelescoPy)
[![Documentation Status](//readthedocs.org/projects/telescopy/badge/?version=latest)](https://telescopy.readthedocs.io/en/latest/?badge=latest)

A telescope observations simulator for Python

## Summary

Telescopy aims to simulate images of astronomical objects observed with optical 
and infrared instruments. It does this by creating models of the optical train 
and astronomical objects and then pushing the object through the optical train. 
The resulting 2D image is then broadcast to a detector chip and read out into a 
FITS file. 

This code was originally based on the [SimCADO](www.univie.ac.at/simcado)

## Dependencies

```
numpy
scipy
astropy
pysynphot
wget
```