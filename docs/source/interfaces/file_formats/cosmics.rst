Cosmic ray file formats
=======================

Cosmic rays maps
----------------

**Description**: A series of images of cosmic ray hits. Only really applicable
to CCD detectors

**File type**: FITS

**File contents**:

* EXT 0 Meta data
* EXT 1 Data in cosmics images
* EXT 2..N Images of Cosmics

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : COSMICS
    EDIM : 2
    ECAT  : 1     # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 2      # In which extension does the real data start


* EXT 1 Header (BinTable)::

    TBD

* EXT 2..N Header (2D image)::

    TBD

**Required data format**

* EXT 1 (BinTable)

  A table containing whatever information is deemed useful to descibe cosmic rays.
  An example might be something like this:

  === ====== ====== =====
  ext energy length angle
  int float  float  float
  ... keV    pixels deg
  === ====== ====== =====

  where:

  * "ext" is the extension number of the image,
  * "energy" is the energy if the cosmic ray that caused the track,
  * "length" is the length of the track on the detector, and
  * "angle" is the rotation angle of the rtack w.r.t to the x-axis.

  .

* EXT 2..N (2D images)

  (x,y) images of various cosmic ray hits.

