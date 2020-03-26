Detector list file format
=========================

**Description**: Describes the physical characteristics of the chips used in
the detector array. The conversion between on-sky coordinates and detector
plane coordinates is handled by the SimCADO parameter ``pixel_scale``

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
    ETYPE : CHIPLIST
    EDIM  : 2

**Required data format**:

An ASCII table with the following columns:

=== ==== ===== ===== ===== ===== ======== ===== ===== ===== ======
id  type x_cen y_cen x_hw  y_hw  pix_size x_len y_len angle gain
--- ---- ----- ----- ----- ----- -------- ----- ----- ----- ------
int str  float float float float float    int   int   float float
... ...  mm    mm    mm    mm    mm       pix   pix   deg   e-/ADU
=== ==== ===== ===== ===== ===== ======== ===== ===== ===== ======

where:

* "id" is a reference id for the chip,
* "type" is the type of chip, e.g. generic_nir / generic_ccd / hawaiiXrg / aquarius,
* "x_cen" and "y_cen" are the physical coordinates of centre of the chip on the
  detector plane in [mm],
* "x_hw", "y_hw" are the half-widths of the chip, i.e. length/2 and height/2 or
  how far the chips extend from the central coordinates,
* "pix_size" is the physical size of pixels in the detector in [mm],
* "x_len", "y_len" are the number pixels in each dimension,
* "angle" is the rotation of the detector relative to the x-axis, and
* "gain" is the conversion factor for electrons (photons) to ADUs


.. WARNING::
    x_len and x_hw are redundant. Discuss which one to keep

