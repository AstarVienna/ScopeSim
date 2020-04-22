Weight map file formats
=======================

Illumination  maps
------------------

**Description**: A single map which can be used to describe the spatial
variation of light over the field of view. The two cases are in spatial
variation of transmission (Illumination map: ILLUMMAP) and emission
(Emission map - EMISMAP).

**File type**: FITS

**File contents**:

* EXT 0 Meta data
* EXT 1 Map

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : ILLUMMAP / EMISMAP
    EDIM : 2
    ECAT  : -1     # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 1      # In which extension does the real data start

* EXT 1 Header (2D image)::

    Standard WCS for the images
    CTYPEn
    CUNITn
    CRVALn   # (0,0) meaning the centre of the field of view
    CRPIXn   # Pixel which corresponds to the centre of the field of view
    CDELTn

**Required data format**

* EXT 1 (2D image)

  An image of the intensity differences over the focal plane. The resolution
  can be much coarser than the detector plate scale. This map will be multiplied
  with a number of photons to represent either the spatial variations in
  background emission, or variation in transmission of a surface / system


Pixel sensitivity maps
----------------------

**Description**: A series of pixel maps for each detector in the instruments
describing the relative sensitivity of each pixel

**File type**: FITS

**File contents**:

* EXT 0 Meta data
* EXT 1 Maps

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : PIXELMAP
    EDIM : 2
    ECAT  : -1     # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 1      # In which extension does the real data start


* EXT 1 Header (2D image)::

    CHIPIDn     # The chip ID for each layer in the data cube, if not sequential

**Required data format**

* EXT 1 (3D image)

  A cube with dimensions (x,y,N) where each (x,y) plane is the pixel sensitivitiy
  map for chip N in the detector array


Persistence maps
----------------

**Description**: A series of maps for each detector in the instrument
describing the persistence image that should be added to each exposure

**File type**: FITS

**File contents**:

* EXT 0 Meta data
* EXT 1 Maps

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : PERSMAP
    EDIM : 2
    ECAT  : -1     # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 1      # In which extension does the real data start


* EXT 1 Header (2D image)::

    CHIPIDn     # The chip ID for each layer in the data cube, if not sequential

**Required data format**

* EXT 1 (3D image)

  A cube with dimensions (x,y,N) where each (x,y) plane is the persistence
  map for chip N in the detector array