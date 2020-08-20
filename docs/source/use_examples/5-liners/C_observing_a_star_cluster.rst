Observing: A star cluster from a template
=========================================
Using the star cluster function from ``scopesim_templates`` with the LFOA telescope

.. jupyter-execute::
    :hide-code:
    :raises:

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    %matplotlib inline

    import scopesim
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = "../temp/"

TL;DR
-----

.. jupyter-execute::
    :raises:

    import scopesim
    import scopesim_templates as sim_tp

    my_cluster = sim_tp.basic.stars.cluster(mass=10000,
                                            distance=2000,
                                            core_radius=1)

    scopesim.server.download_package("telescopes/LFOA")
    lfoa = scopesim.OpticalTrain("LFOA")
    lfoa.observe(my_cluster)
    lfoa.readout(filename="TEST.fits")


This 5-liner uses concepts from:

- :doc:`A_loading_packages`: downloading instrument packages

See also:

- :doc:`B_source_from_arrays`: downloading instrument packages
- :doc:`B_source_from_image`: downloading instrument packages
- :doc:`B_source_from_table`: downloading instrument packages

Explanation
-----------

The star cluster templates function
+++++++++++++++++++++++++++++++++++
`ScopeSim Templates <https://scopesim-templates.readthedocs.io/en/latest/>`_ is a support package for ``ScopeSim``.
Currently this package only contains a few basic helper functions for creating standard on-sky targets in the ``scopesim.Source`` format.
One of these basic functions creates a star cluster from a set of physical parameters:

.. jupyter-execute::
    :raises:

    my_cluster = sim_tp.basic.stars.cluster(mass=10000, distance=2000, core_radius=1)

Here we create a cluster with 10000 solar masses, at a distance of 2000 parsec, with a core radius of 1 parsec.
The star masses are drawn from a Kroupa (2001) IMF, the positions are drawn from a Gaussian distribution, and the associated spectra for each star is from the Pickles (1998) catalogue.

The table describing the cluster is held in the ``src.fields`` list:

.. jupyter-execute::
    :raises:

    print(my_cluster.fields[0])

Here we can see that the ``scopesim_templates`` function has done all the hard work of working out the on-sky size, distance modulus, and geometric parameters of the cluster.

.. note:: Community code contributions to ``scopesim_templates`` are most welcome!

    If you have code that creates a spectro-spatial description of an on-sky object, and you would like this object to be included in ``scopesim_templates``, please make a pull request via the `ScopeSim_Templates GitHub repository <https://github.com/astronomyk/scopesim_templates/pulls>`_


The Leopold-Figl Observatory for Astrophysics
+++++++++++++++++++++++++++++++++++++++++++++

The LFOA is the 1.5m telescope that belongs to the `Department of Astrophysics at the University of Vienna <https://foa.univie.ac.at/>`_.
The telescope's camera has 1092 x 736 pixels, covering a 5.58 x 3.75 arcminute field of view.

We download the LFOA package using the standard method from :doc:`A_loading_packages`:

.. jupyter-execute::
    :raises:

    scopesim.server.download_package("telescopes/LFOA")

To simply observe using default telescope values, we can use the shortcut option and create an optical model directly:

.. jupyter-execute::
    :raises:

    lfoa = scopesim.OpticalTrain("LFOA")

If we want to set more andvaced features, like selecting a different filter, we need create a ``UserCommands`` object, and set the bang-string keyword ``!OBS.filter_name``:

.. jupyter-execute::
    :raises:

    cmds = scopesim.UserCommands(use_instrument="LFOA")
    cmds["!OBS.filter_name"] = "sloan_z"
    lfoa = scopesim.OpticalTrain(cmds)

As a side note, if the sky background is too low, we can also increase this with the bang-string keyword ``!OBS.sky.bg_mag``.

.. note:: Top-level control parameters are contained in a ``UserCommands`` object.

    If we have an external ``UserCommands`` object, these can be viewed by simply printing the objects::

        print(cmds)

    If we have already built an optical model, these commands are contained in ``<OpticalTrain>.cmds``.
    For the LFOA these can be viewed by calling ``print(lfao.cmds)``

We can view the spectral response of the system by using internal optic manager:

.. jupyter-execute::
    :raises:

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    %matplotlib inline

    wave = np.arange(3000, 11000)      # in Angstrom, [default units of SynPhot]
    thru = lfoa.optics_manager.surfaces_table.throughput(wave)

    plt.plot(wave, thru)
    plt.show()


Observing the cluster with the Figl telescope
+++++++++++++++++++++++++++++++++++++++++++++

To observe the cluster with the LFOA telescope, we call the ``observe`` method and pass the source object:

.. jupyter-execute::
    :raises:

    lfoa.observe(my_cluster)

This generates an "expectation" image on the image plane directly above the detector in units of ``ph/s/pixel``.
This image contains no noise.
It is used as the basis for generating the detector readout image.

.. jupyter-execute::
    :raises:

    im = lfoa.image_planes[0].image
    plt.imshow(im, norm=LogNorm())

To make the "raw" data for the telescope, we call the ``readout`` method.
We can provide a ``filename`` if we want to save a ``FITS`` image to disc;

.. jupyter-execute::
    :raises:

    lfoa.readout(filename="TEST.fits")

Or we can work directly with the returned list of ``astropy.fits.HDUList`` objects:

.. jupyter-execute::
    :raises:

    hdus = lfoa.readout()

Here we must be careful though. ``ScopeSim`` returns a list of ``HDUList`` objects, not just a single one, even though there is only one detector on the Figl observatory.
This is because the software is set up to simulate instruments with multiple detector arrays (e.g. XSHOOTER). To avoid differing API endpoints for different instruments, the decision was made to always return a list, even if there is only one detector in the instrument.

.. jupyter-execute::
    :raises:

    im = hdus[0][1].data
    plt.imshow(im, norm=LogNorm())
    plt.colorbar()


Updating the exposure time
++++++++++++++++++++++++++

The exposure time (``dit``, and/or ``ndit``) are dynamical parameters and do not require the optical model to be remade.
Hence these can be updated at any point using the ``.cmds`` command object inside the telescope model:

.. jupyter-execute::
    :raises:

    lfoa.cmds["!OBS.dit"] = 1
    hdus = lfoa.readout()

    im = hdus[0][1].data
    plt.imshow(im, norm=LogNorm())
    plt.colorbar()
