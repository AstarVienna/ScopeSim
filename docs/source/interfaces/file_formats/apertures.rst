Aperture list file formats
==========================

Spectroscopic aperture list
---------------------------

**Description**: Describes the spatial on-sky characteristics for spectrographic
apertures. E.g. which part of the sky the fibres of a MOS see, or which parts
of the sky the pseudo-slits of an image-slicer IFU see.

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
    ETYPE : APERLIST
    EDIM  : 2

**Required data format**:

.. warning:: Out of date!

An ASCII table with the following columns:

=== ==== ====== ====== ====== ====== ===== ======
id  type dra    ddec   hw1    hw2    angle s_off
--- ---- ------ ------ ------ ------ ----- ------
int str  float  float  float  float  float float
... ...  arcsec arcsec arcsec arcsec deg   arcsec
=== ==== ====== ====== ====== ====== ===== ======

where:

* "id" is the number of the aperture,
* "type" is slit / fibre,
* "dra", "ddec" are the position of the aperture relative to the centre of
  the field of view in arcsec,
* "hw1", "hw2" are the half-widths of the aperture in arcsec.
  For a slit aperture these refer to half the length and half the width
  (e.g. a 15" x 1" slit would have hw1=7.5" and hw2=0.5").
  For a fibre aperture these refer to the radii of the semi-major and semi-minor
  axes. If the aperture is perfectuly circular then hw1==hw2.
* "angle" is the angle of rotation of the slit or fibre w.r.t to the RA axis,
* "s_off" (relevent only for slits) is the positional offset along the slit of
  the refenence position. E.g. if the trace desciption is not symetrical and
  requires an offset.