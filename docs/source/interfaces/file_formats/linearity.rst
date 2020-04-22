Detector Linerarity file formats
================================

Linearity Curve
---------------

**Description**: The relationship between actual photon counts and registered
photon counts for a detector chip

**File type**: ASCII

**File contents**:

* Header info, commented out with either "#" or "\"
* ASCII table

**Required header keywords**::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : LINEARIT
    EDIM  : 0

**Required data format**:

An ASCII table with the following columns:

========= ===============
real_flux detected_flux
--------- ---------------
int       int
photons   photo-electrons
========= ===============

where:

* "real_flux" is the real incoming photon flux per pixel
* "detected_flux" is the photon flux per pixel registered by the detector chip
