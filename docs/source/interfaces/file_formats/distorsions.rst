Distortion map file formats
===========================

2D Discretised distortion maps
------------------------------

**Description**: Maps which cover the extent of an image plane and describe the
extent of the distortion in both spatial dimensions. The structure allows for
wavelength dependent distortions.

**File type**: FITS

**File contents**:

* EXT 0 Meta data
* EXT 1 Catalogue
* EXT 2..N Distortion maps (x,y,2)

**Required header keywords**:

* EXT 0 Header (Empty)::

    AUTHOR
    DATE
    ORIGDATE
    SOURCE
    STATUS
    ETYPE : DISTMAP
    EDIM : 5
    ECAT  : -1     # In which extension is the catalogue data. -1 if no catalogue
    EDATA : 1      # In which extension does the real data start

* EXT 1..N Header (3D image)::

    WAVE0       # Wavelength valid for extension. -1 if achromatic
    WAVEUNIT    # Unit of wavelength. If absent assumption is [um]
    Standard WCS for the images
    CTYPEn
    CUNITn
    CRVALn   # (0,0) meaning the centre of the field of view
    CRPIXn   # Pixel which corresponds to the centre of the field of view
    CDELTn

**Required data format**

* EXT 1..N (3D image)

  N-1 Image cubes (x,y,2)

  The 2 layers of the cube will describe the amount of distortion in each of the
  x and y dimensions over the field. If there is wavelength dependent distortion,
  each extension describes the distortion valid for the wavelength definied by
  the WAVE0 keyword in the header.

