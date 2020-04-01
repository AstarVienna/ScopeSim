Surface list file format
========================

**Description**: A list of surfaces and links to the files which contain the
spectral characteristics for Transmission, Emission, and Reflection.

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
    ETYPE : SURFLIST
    EDIM  : 1

**Required data format**:

An ASCII table with the following columns:

===== ==== ==== ===== ===== ===== =========== ========
order name type outer inner angle temperature filename
----- ---- ---- ----- ----- ----- ----------- --------
int   str  str  float float float float       str
...   ...  ...  m     m     deg   degC        ...
===== ==== ==== ===== ===== ===== =========== ========

where:

* "order" is the position along the optical path, i.e. M1 is 1, M5 is 5,
* "name" of the element,
* "type" of surface regarding throughput: reflective (r) or transmittive (t)
* "outer", "inner" are the outer and inner diameters in [m] of the optical element,
* "angle" is the angle at which the element is rotated w.r.t the optical axis,
* "temperature" is the temperature in degrees Celcius of the optical element,
* "filename" refenences the file containing the spectral response curves for
  transmission, emission, and reflection
