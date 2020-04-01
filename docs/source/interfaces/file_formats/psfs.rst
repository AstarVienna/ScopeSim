PSFs file formats
=================

Field varying PSFs
------------------
**Description:** Contains PSFs for 1..N wavelengths that are applicable for a
small part of the field. Each extension should contain a cube with each layer
containing the PSF for a different region in the field of view. The applicable
on-sky region can either be described in EXT 1 in table format, or with a
weights map that covers the field of view

**File type**: FITS

**File contents**:

* Ext 0 contains meta data,
* Ext 1 contains either a table or a weights map
* Ext 2..N contains the wavelength dependent PSF cubes

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : FVPSF
    EDIM  : 5
    ECAT  : 1      # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 2      # In which extension does the real data start

* EXT 1 Header (BinTable)::

    NUMPSFS : 1..M   # How many PSF layers per cube.
    CATTYPE : table  # Catalogue format used to describe the valid FOV for a PSF
    CUNIT1           # Units for values in table (arcsec / arcmin / deg)

* EXT 1 Header (Image 2D/3D)::

    NUMPSFS : 1..M   # How many PSF layers per cube.
    CATTYPE : image
    Standard WCS for the image
    CTYPEn
    CUNITn
    CRVALn   # (0,0) meaning the centre of the field of view
    CRPIXn   # Pixel which corresponds to the centre of the field of view
    CDELTn

* EXT 2..N Header (Image 3D)::

    WAVE0
    WAVEUNIT    # Unit of wavelength. If absent assumption is [um]
    Standard WCS for the image
    CTYPEn
    CUNITn
    CRVALn   # (0,0) meaning the centre of the field of view
    CRPIXn   # Pixel which corresponds to the centre of the field of view
    CDELTn

**Required data format**

* EXT 1 (BinTable)

  (N,3) Table with the following columns

  ====== ====== =====
  x      y      layer
  float  float  int
  arcsec arcsec none
  ====== ====== =====

  where:

  * "x","y" are the centres of the valid regions. SimCADO draws its own map to
    define where the borders are between these regions
  * "layer" is the position along the M dimension of the PSF cube

  .

* EXT 1 (Image 2D/3D)

  Image cube (x,y,N) with N>=1 layers

  Each layer is an image of the whole focal plane (can use much coarser
  resolution than plate scale) where the pixel values correspond to the PSF layer
  (from EXT >=2) that should be used in a given region. If there is a different
  map for each wavelength then number of layers in the EXT 1 cube should equal to
  the number PSF extensions, i.e. size(EXT 1) = (x,y,N) for a file with N
  extensions. If there is only one layer in EXT 1, it will be assumed that this
  weight map works for all wavelengths.

* EXT 2..N (Image 3D)

  N-2 Image cubes (x,y,M) each with M>=1 layers

  Each EXT holds a cube with PSFs for a certain wavelength. Each layer (x,y) in a
  cube is a PSF kernel which is valid for a certain region of the focal plane and
  for the wavelength given by WAVE0 in the header of each EXT. The location of
  the valid region is given by the data in EXT 1.


Spatially constant PSFs
-----------------------
**Description**: Contains PSFs for 1..N wavelengths that are applicable over
the whole field of view. The data structure will be the same as the field
varying PSFs.

**File type**: FITS

**File contents**:

* Ext 0 contains meta data,
* Ext 1 contains either a table or a weights map
* Ext 2..N contains the wavelength dependent PSF cubes

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : CONSTPSF
    EDIM : 3
    ECAT  : 1      # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 2      # In which extension does the real data start

* EXT 1 Header (Empty)::

    NUMPSFS : 1
    CATTYPE : none  # The type of data used to describe the valid FOV for a PSF

* EXT 2..N Header (Image 2D)::

    WAVE0
    WAVEUNIT    # Unit of wavelength. If absent assumption is [um]
    Standard WCS for the image
    CTYPEn
    CUNITn
    CRVALn   # (0,0) meaning the centre of the field of view
    CRPIXn   # Pixel which corresponds to the centre of the field of view
    CDELTn

**Required data format**

* EXT 1 (Empty)

  No data unit needs to be attached

* EXT 2..N (Image 2D/3D)

  N-2 Image cubes (x,y,1)

  Each EXT holds a PSF kernel for a certain wavelength given by WAVE0 in the
  header of each EXT.
