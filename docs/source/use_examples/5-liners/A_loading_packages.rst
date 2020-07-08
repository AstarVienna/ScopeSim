Control: Loading Packages
=========================

TL;DR
-----

.. jupyter-execute::
    :raises:

    import scopesim, os
    if not os.path.exists("../temp/"): os.mkdir("../temp/")

    scopesim.rc.__config__["!SIM.file.local_packages_path"] = "../temp/"

    pkg_names = ["locations/Paranal", "telescopes/VLT", "instruments/HAWKI"]
    scopesim.server.download_package(pkg_names)

    cmds = scopesim.UserCommands(use_instrument="HAWKI")
    hawki = scopesim.OpticalTrain(cmds)


Explanation
-----------

Before we can load anything we need to download the instrument packages from the
instrument reference database (https://github.com/astronomyk/irdb).

We can list all the available packages like so:

.. jupyter-execute::
    :raises:

    scopesim.server.list_packages()


Note that the packages are split into three categories:

* locations
* telescopes
* instruments

To use an instrument package, we will need to download the support packages (location and telescope) that are relevant to the instrument.
In the case HAWKI, this means also getting the VLT and Paranal packages:

.. jupyter-execute::
    :raises:

    pkg_names = ["locations/Paranal", "telescopes/VLT", "instruments/HAWKI"]
    scopesim.server.download_package(pkg_names)


The standard way to load an instrument package is to create a ``UserCommands``
object and the ``use_instrument=`` parameter:


.. jupyter-execute::
    :raises:

    hawki_cmds = scopesim.UserCommands(use_instrument="HAWKI")

This will allow us to play around with the parameters before we generate a model
of the whole optical system.
The optical model is created by passing the ``UserCommands`` object to an
``OpticalTrain`` object:

.. jupyter-execute::
    :raises:

    hawki = scopesim.OpticalTrain(hawki_cmds)

However if we are happy to accept all the default values and simply want to
simulate an observation, we can bypass the user commands step, and initialise
the ``OpticalTrain`` with the name of package that we have on the local disc.

.. jupyter-execute::
    :raises:

    hawki = scopesim.OpticalTrain("HAWKI")

.. note:: Packages are stored by default in the working directory

   This can be changed by setting the following ``rc`` entry::

       scopesim.rc.__config__["!SIM.file.local_packages_path"]

