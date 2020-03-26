Transmission, Emissivity, Reflection (TER) Curve file format
============================================================

**Description**: A table containing the wavelength dependent coefficients for
the spectral response of an optical element.

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
    ETYPE : TERCURVE
    EDIM  : 1

**Optional keywords**::

    MATERIAL     # For optical element coatings or substrates
    AIRMASS      # For atmospheric TER curves
    PWV          # For atmospheric TER curves

**Required data format**:

An ASCII table with the following columns:

===== ============ ========== ==========
lam   transmission emissivity reflection
----- ------------ ---------- ----------
float float          float      float
um    0..1           0..1       0..1
===== ============ ========== ==========

where

* "lam" is the wavelength in [um],
* "transmission" is the coefficient of transmission between [0,1]
* "emissivity" is the coefficient of emissivity between [0,1]
* "reflection" is the coefficient of reflection between [0,1]

In general the transmission + reflection should equal 1. Emissivity is a
the coefficient applied to a blackbody emission curve for the optical element.
