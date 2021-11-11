---
title: 'ScopeSim - A pythonic astronomical instrumental data simulation engine'
tags:
  - Python
  - Astronomy
  - Simulations
  - Telescopes
  - Instruments
  - Extreme Large Telescope

authors:
  - name: Kieran Leschinski
    orcid: 0000-0003-0441-9784
    affiliation: 1
  - name: Oliver Czoske
    orcid: 0000-0003-3127-5341 
    affiliation: 1
  - name: Miguel Verdugo
    orcid: 0000-0001-5027-557X
    affiliation: 1 
  - name: Hugo Buddelmeijer
    orcid: 0000-0001-8001-0089
    affiliation: 2
  - name: Gijs Verdoes-Kleijn
    orcid: 0000-0001-5803-2580
    affiliation: 2
  - name: Werner Zeilinger
    orcid: 0000-0001-8507-1403
    affiliation: 1
  - name: Joao Alves
    orcid: 0000-0002-4355-0921
    affiliation: 1
    
affiliations:
 - name: Department of Astrophysics, University of Vienna
   index: 1
 - name: OmegaCEN, Kapteyn Astronomical Institute, University of Groningen
   index: 2
   
date: 28 September 2021
bibliography: paper.bib

---

# Summary

- A pythonic simulation engine for astronomical instrument data products
- It 


Documentation can be found at https://scopesim.readthedocs.io/en/latest/

# Statement of need

- Why we need ScopeSim
    - Each consortium invests time and effort in writing simulators specifically for their instrument
    - Once the commisioning of the instrument is done, the simulator is forgotten
    - At any one time there are few instruments being built, thus no effort has gone into keeping code and knowledge    
    - The majority of astronomical instruments contain the same optical elements
    - There is no standard interface for desribing instrumental effects and no standard code library (like astropy) 
    - The ScopeSim framework provides the building blocks that each simulator needs, thus eliminating the need to start from scratch
    - With a standard simulation engine for multiple instruments, it becomes much easier to make meaningful comparisons between output data. Compare apples to apples

- Audiences
    - Scientists, feasibility studies
    - Scientists, observation proposals  
    - Data redcution pipeline developers
    - New PIs, Proposals for new instruments

![caption](path)

# ScopeSim workflow

## Connection to other packages in the software framework

## Basic code example 





# Acknowledgments

ScopeSim depends on the following packages: 
Numpy [@numpy],
SciPy
Astropy [@astropy2018].
SynPhot

This project was funded by project IS538004 of the Hochschulraum-strukturmittel (HRSM) provided by the Austrian Government and administered by the University of Vienna.

# References
