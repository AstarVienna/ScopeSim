Spectral Traces File Format
===========================

**Description:** A file to hold all of the trace maps for a spectrograph.
The catalogue (EXT 1) connects a single trace to a sky mask (either slit or
fibre) and a specific image plane. The traces provide the position on the image
plane where light of a certain wavelength will fall.

.. Caution::
    Dependency Warning

    This number of masks references in the EXT 1 table must be compatible with
    the number of masks described in the file describing masks. These are
    separate files as the positions of the masks (in the case of fibres) is not
    always fixed on the focal plane. Hence different spatial configurations
    for MOS instruments should still reference the same trace layout.

**File type**: FITS

**File contents**:

* EXT 0: Meta data
* EXT 1: BinTable with the catalogue connecting projection traces to on-sky
  apertures and image planes
* EXT 2..N: BinTables, each with a single trace mapping in the detector plane

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : SLITTRAC / FIBRTRAC
    EDIM : 3
    ECAT  : 1       # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 2       # In which extension does the real data start

* EXT 1 Header (BinTable)::

    CATTYPE : table
    NUMMASKS : 1..M

* EXT 2..N Header (BinTable)::

    Optional header info
    TRACNAME    # Name of trace
    MASKID      # Mask number that produces the trace
    PLANEID     # Image plane number - which image plane receives the trace

**Required data format**

* EXT 1 (BinTable)

  The catalogue table which connects a trace to a mask (slit/fibre). It should
  contain the following columns:

  === ==== =========== =============
  ext name aperture_id imageplane_id
  --- ---- ----------- -------------
  int str  int         int
  === ==== =========== =============

  where:

  * "ext" is the extension number (2..N) for the Trace,
  * "name" is the name of the fibre / order / slit,
  * "aperture_id" is identifying number of the slit / fibre in the file
    containing the description of the slits / fibres
  * "imageplane_id" is the identifying number of the image plane that will be
    used for the projection of the aperture


* EXT 2 (BinTable)

  Each extension represents the path of a single spectral trace over the detector.
  There are two types of trace: SLITTRAC / FIBRTRAC.
  A slittrace preserves the spatial extent of the incoming light and thus
  describes a series of lines or curves that the projected slit mask will follow
  for different wavelengths. The number of points per wavelength used to trace
  the projection of the slit is unlimited (in theory). SimCADO will recognise the
  width of the table and calculate how many points are included for the
  polynomial fit

  For a SLITTRAC, the table should contain the following columns:

  ====== ====== ====== ====== === ====== ====== ======
  wave   s1     x1     y1     ... sN     xN     yN
  ------ ------ ------ ------ --- ------ ------ ------
  float  float  float  float  ... float  float  float
  micron arcsec mm     mm     ... arcsec mm     mm
  ====== ====== ====== ====== === ====== ====== ======

  where:

  * "lam" is wavelength,
  * "s" is the position along the slit relative to the reference point of the mask
    (defined in the mask description file),
  * "x", "y" are the positions on the detector plane (in mm) of each wavelength.

  For a FIBRTRAC, there is no slit dimension, as the fibre scrambles the spatial
  structure of the incoming light. Instead the exiting beam has a certain width.
  The trace table should contain the following columns:

  ====== ====== ====== ====== ====== ======
  lam    x      y      dx     dy     ang
  ------ ------ ------ ------ ------ ------
  float  float  float  float  float  float
  micron mm     mm     mm     mm     deg
  ====== ====== ====== ====== ====== ======

  where:

  * "lam" is wavelength,
  * "x", "y" are the positions on the detector plane (in mm) of each wavelength,
  * "dx", "dy" are the width and height of the fibre beam projected on the focal
    plane, and
  * "ang" is the rotation angle of the projected fibre beam w.r.t to the x axis
