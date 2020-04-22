Spectral Response Curve file format
===================================

Throughput Curves
-----------------
E.G. Filter Curve, TER Curve (Transmission, Emissivity, Reflection)

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


Emission curves
---------------

**Description**:

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
    ETYPE : EMSCURVE
    EDIM  : 1
    EUNIT       # units of emission e.g. ph s-1 m-2 arcsec-2 bin-1

**Required data format**:

An ASCII table with the following columns:

===== ========
lam   emission
----- --------
float float
um    EUNIT
===== ========

where:

* "lam" is the wavelength in [um]
* "emission" is the wavelength dependent emission at the given wavelength. The
  units of the emission are defined by EUNIT in the header