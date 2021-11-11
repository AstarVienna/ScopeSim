---
title: 'AnisoCADO: a python package for analytically generating adaptive optics point spread functions for the Extremely Large Telescope'
tags:
  - Python
  - astronomy
  - simulations
  - point spread functions
  - Extreme Large Telescope
authors:
  - name: Kieran Leschinski
    orcid: 0000-0003-0441-9784
    affiliation: 1
  - name: Eric Gendron
    affiliation: 2 
affiliations:
 - name: Department of Astrophysics, University of Vienna
   index: 1
 - name: Observatoire de Paris
   index: 2
date: 15 June 2020
bibliography: paper.bib
---

# Summary

AnisoCADO is a Python package for generating images of the point spread function (PSF) for the european extremely large telescope (ELT). 
The code allows the user to set a large range of the most important atmospheric and observational parameters that influence the shape and strehl ratio of the resulting PSF, including but not limited to: the atmospheric turbulence profile, the guide star position for a single conjugate adaptive optics (SCAO) solution, differential telescope pupil transmission, etc.
Documentation can be found at https://anisocado.readthedocs.io/en/latest/


# Statement of need

## Adaptive optics are mandatory for the next generation of ground-based telescopes
The larger the telescope aperture, the smaller the diffraction limit of the observations. 
For space-based telescope this statement is always true. 
However the resolution of ground based telescopes is limited by the blur caused by turbulence in the atmosphere - known as atmospheric Seeing. 
This blurring can be (mostly) removed by measuring the deformation of the wavefront of the incoming light, and applying an equal and opposite deformation to the surface of one or more of the mirrors along a telescope's optical path.
The current fleet of large (8-10m) telescopes were built to primarily operate at the edge of the natural seeing limit (FWHM~0.5 arcseconds @ 1um). 
Over the last two decades some have received upgrades in the form of active and adaptive mirrors in order to achieve up to 20x increase in resolution afforded by the physical diffraction limit of a 10m primary mirror (FWHM~0.03 arcseconds @ 1um).
The next generation of "extremely large" telescopes will have primary mirrors on the order of 30-40m, with theoretical diffraciton limits on the order of 50x smaller than the natural Seeing limit.
In order for these telescopes to resolve structures at scales of the diffraction limit, they must, by design, include adaptive optics systems.

## Diffraction limited point-spread-functions are complex beasts
The point spread function (PSF) of an optical system is the description of the spatial distribution of light from an infinitely small point source after passing through an optical system (e.g. layers of the atmosphere, mirrors of a telescope).
Due to the random nature of atmospheric turbulence, the PSF of a star in a Seeing-limited observation is well approximated by a ("nice") smooth Gaussian-like function.
The PSF of a diffraction limited telescope system using an adaptive-optics correction is a complex ("ugly") function that depends on a veritible zoo of atmospheric, observational, and technical parameters. 
From an astronomers point of view, the consequences of a poor adaptive optics solution means the difference between a successful and a failed observation run.
Therfore it is imperative that the consequences of such large variations in the PSF are accounted for in advance by those proposing to observe with this next generation of billion-dollar telescopes.


## AnisoCADO - Anisoplanatism for MICADO

AnisoCADO (Anisoplanatism for MICADO) is a package for generating images of the point spread function for a given set of atmospheric, observational, and technical conditions.
It does this by combining a series of wavefront phase screens from the elements of the atmosphere and AO system that influence the final AO correction (e.g. atmospheric anisoplantism, WFS aliasing, actuator fit, etc). 
The final phase screen is applied to the optical transfer function for the telescope optical system. 
The resulting image is the expected PSF for a long exposure (>10s) on-axis observation at the given wavelength.
For single-conjugate adaptive optics modes, the field PSF degrades as distance from the guide star increases. 
This effect is taken into account by shifting the anisoplanatic phase screen relative to the calculated phase screen correction for the deformable mirror. 
Figure \autoref{fig:psf_grid} shows how the PSF changes with distance from an on-axis guide star.   
For a more detailed discussion of the mathematics behind anisoplanatism in the context of the ELT, the reader is referred to @clenet2015.

![A grid of Ks-band (2.15um) PSFs for a range of distances from the natural guide star. 
The PSFs were generated using the ESO median turbulence profile. 
\label{fig:psf_grid}](Ks-band_psf_grid.png)


### Inputs
The final ELT PSF is the combination of many factors. The vast majority of these are irrelevant for the casual user. 
AnisoCADO therefore provides three preset option, corresponding to the standard ESO Q1, Median and Q4 turbulence profiles.
All other parameters are initialised with default values.
For the case of a SCAO system (for which AnisoCADO was originally conceived) PSFs can be generated for multiple guide star offsets without needing to re-make all phase screens by using the special class method ``.shift_off_axis(dx, dy)``  

For more detailed use cases, the following parameters are available to the user:

| Atmosphere                    | Observation                      | Telescope                    |
|-------------------------------|----------------------------------|------------------------------|
| * turbulence profile          | * natural guide star position    | * pupil image                |
| * height of turbulent layers  | * central wavelength             | * 2D pupil transmissivity    |
| * stregth of turbulent layers | * pupil rotation angle           | * dead/empty mirror segments |
| * wind speed                  | * Zenith distance of observation | * plate scale                |
| * Seeing FWHM @ 500nm         |                                  | * residual wavefront errors  |
| * Fried parameter             |                                  | * AO sampling frequency      |
| * outer scale                 |                                  | * AO loop delay              |
|                               |                                  | * Interactuator distance     |



### Outputs

AnisoCADO is easily integrated into the standard astronomers toolbox. 
PSF images generated by AnisoCADO can be output as either ``numpy`` arrays, or standard ``astropy.io.fits.ImageHDU`` objects. 
The latter can be written to file using the standard ``astropy`` syntax.

As AnisoCADO was written to support the development of the MICADO instrument simulator [@simcado2016; @simcado2019], it is also possible to generate ``FieldVaryingPSF`` objects using the helper functions in the ``misc`` module.
Such files are also compatible with the generic instrument data simulator framework, ScopeSim [@scopesim].


Basic Example
-------------
The AnisoCADO API is described in the online documentation, which can be found at: <https://anisocado.readthedocs.io/>. For the purpose of illustration, the following 5 lines were used to generate the grid of PSFs in figure \autoref{fig:psf_grid}.

```
import numpy as np
from anisocado import AnalyticalScaoPsf

psf = AnalyticalScaoPsf()
psf_grid = []
for x, y in np.mgrid[-14:15:7, -14:15:7].flatten().reshape((2, 25)).T:
    psf.shift_off_axis(x, y)
    psf_grid += [psf.kernel]
```


# Acknowledgments

AnisoCADO depends on the following packages: 
Numpy [@numpy],
Matplotlib [@numpy],
Astropy [@astropy2018].

This development of this project was funded by the project IS538004 of the Hochschulraum-strukturmittel (HRSM) provided by the Austrian Government and administered by the University of Vienna.


# References


@misc{scopesim,
    author = {{Leschinski}, Kieran},
     title = "{ScopeSim - A python framework for creating astronomical instrument data simulators}",
      year = {2020},
 publisher = {​GitHub},
   journal = {​GitHub repository},
       url = {​https://github.com/astronomyk/scopesim}
}