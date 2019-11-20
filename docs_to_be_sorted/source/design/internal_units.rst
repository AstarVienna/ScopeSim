Internal Units
==============
In order to avoid the slow down incurred by using astropy units, ScopeSim uses
the following internal units:

Spectra
-------
* wavelength: ``um`` - when invoking a ``Synphot`` routine, this is converted
  to ``Angstrom``
* emission: ``PHOTLAM`` (``ph s-1 AA-1 cm-2``) - to be consistent with ``Synphot``
* surface emission: ``PHOTLAM arcsec-2``
* transmission, reflection, emissivity: ``unitless`` values in the range [0..1]

On-sky distances
----------------
* x, y, left, right, top, bottom: ``arcsec`` from optical axis
* angle: ``degrees`` relative to zenith angle
* plate_scale: ``arcsec mm-1``
* pixel_scale: ``arcsec pixel-1``

Detector-plane distances
------------------------
* x, y, xhx, yhw: ``mm``
* pixsize: ``mm``
* gain: ``ph ADU-1``
* dark current: ``electron s-1``
* read noise: ``electron``

Mirror characteristics
----------------------
* outer: ``m``
* inner: ``m``
* angle: ``degrees`` relative to the optical axis
* temperature: ``deg_C``
* wavefront_error: ``nm``








