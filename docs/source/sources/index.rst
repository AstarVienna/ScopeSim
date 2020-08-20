Source (target) interface
=========================

.. note:: This page described the interface for the ``Source`` object.

    For information on the helper function which create ``Source`` objects for specific astronomical targets (e.g. galaxy, star cluster), please see the `ScopeSim Templates`_ documentation.


A ScopeSim ``Source`` object is a 2+1D (x, y, lambda) description of an on-sky object.
Hence to build one of these objects we need the following data:

* A spatial description

  Either in table form (point sources) or in image/bitmap form (extended source)

* A spectral description

  Basically a spectrum with wavelength and flux information

For example we could use two FITS images to describe the spatial flux distribution of the young and old components of a spiral galaxy, and two spectra, e.g. from a starburst and an elliptical galaxy, for the spectral description.

Internally this information is stored in the following attributes:

* ``<Source>.fields`` - the spatial information [Table or ImageHDU from astropy]
* ``<Source>.spectra`` - the spectral information [SourceSpectrum from synphot]

More detailed information about the ``scopesim.Source`` object can be found in the ``scopesim`` documentation.

Below is a description of how to initialise (create) a ``Source`` object.


Spatial description
-------------------

Series of arrays
++++++++++++++++
The most basic method for creating point source ``Source`` objects is by
passing a series of arrays to the following argument keys:

* ``x, y`` : positions [arcsec] from the centre of the field of view
* ``ref`` : refers to the spectrum associated with a given point source

Optionally, if multiple point sources use the same spectrum, a scaling factor
can be applied to each of the degenerate point sources:

* ``weight`` : scaling factor

Code Example::

   src = Source(x=[1, 2], y=[-1, 5], ref=[0, 0], weight=[3.14, 0.5])

Astropy Table
+++++++++++++


FITS ImageHDU
+++++++++++++
Extended souces should provide an intensity map in the form of an
``fits.HDUImage`` object. The header keywords should contain information
regarding position or the centre of the image (``CRPIXn``) relative to the
centre of the field of view (``CRVALn``) [degrees], and pixel scale
[``CDELTn`` | ``CDn_m``] [degress / pixel].

Each ``ImageHDU`` should also be
associated with an appropriate spectrum (``synphot.SourceSpectrum``) contained
in the list ``<Source>.spectra``. The keyword ``SPEC_REF`` must be present and
refer to the list index of the corresponding spectrum in ``<Source>.spectra``.


Spectral description
--------------------




Series of arrays
++++++++++++++++
* ``wave``
* ``flux``


``synphot.SourceSpectrum``
++++++++++++++++++++++++++



.. _ScopeSim:    https://scopesim.readthedocs.io/en/latest/
.. _`ScopeSim Templates`: https://scopesim-templates.readthedocs.io/en/latest/
.. _IRDB:        https://github.com/astronomyk/irdb
.. _AnisoCADO:   https://anisocado.readthedocs.io/en/latest/
.. _skycalc_ipy: https://skycalc-ipy.readthedocs.io/en/latest/
.. _Pyckles:     https://scopesim-templates.readthedocs.io/en/latest/