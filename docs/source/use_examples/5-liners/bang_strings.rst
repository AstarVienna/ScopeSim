Using Bang (!) strings to control ScopeSim
==========================================

TL;DR
-----

.. jupyter-execute::
    :hide-code:

    import os, scopesim
    pkg_path = os.path.join(os.getcwd(), "temp")
    if not os.path.exists(pkg_path):
        os.mkdir(pkg_path)
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = pkg_path

.. jupyter-execute::

    hawki = scopesim.OpticalTrain("HAWKI")

    hawki.cmds["!ATMO"]
    hawki.cmds["!SIM.random.seed"] = 9001

    hawki.effects
    hawki["filter_curve"].meta["filter_name"]
    hawki.cmds["!OBS.filter_name"] = "H"


This 5-liner uses concepts from:
- :doc:`loading_packages`: downloading instrument packages


Explanation
-----------

Top level parameters
++++++++++++++++++++

Let's assume we are in out working directory and have already downloaded the
packages needed to model HAWKI at the VLT (i.e. Paranal, VLT, HAWKI).
If not, see :doc:`loading_packages`

We can start by loading an ``OpticalTrain`` object for HAWKI

.. jupyter-execute::

    hawki = scopesim.OpticalTrain("HAWKI")

The commands are located in ``hawki.cmds`` and the ``Effect``-objects are
located deep within the depths of the ``OpticalTrain``. They can however be
retrieved if you know what you are looking for.

To view all "user"-facing commands, we can simply print ``hawki.cmds``.
To print just a subset, use "!" plus the alias. E.g:

.. jupyter-execute::

    hawki.cmds["!OBS"]

.. note:: Aliases are quick ways to refer to command catagories

    As ScopeSim treats independent parts of the optical train separately, they
    can each be given a different alias (set in the ``yaml`` files).
    The standard aliases are (with obvious meanings):

    ``!ATMO``, ``!TEL``, ``!RO`` (relay optics), ``!INST``, ``!DET``

    Additionally there are two special categoies ``!OBS`` and ``!SIM`` for the
    observation and simulation parameters respectively.

As ScopeSim uses a YAML structure (nested dictionaries) for the configuration
files, some parameters may themselves also be dictionaries (of dictionaries).
We can navigate down the layers using the "." separator:

.. jupyter-execute::

    hawki.cmds["!SIM.random.seed"] = 9001


Lower level parameters
++++++++++++++++++++++

The top level parameters should contain all the levers the casual user may want
to play with.
If we want to control effects that would normally be hidden from an observer, we
need to know the name of the effect we are looking for.

To list all the effects contained in the HAWKI system, we call:

.. jupyter-execute::

    hawki.effects

By treating ``hawki`` as a dictionary, we can access the individual ``Effect``
objects. The configuration parameters are contained in the ``.meta`` dictionary.

.. jupyter-execute::

    hawki["filter_curve"].meta["filter_name"]

Here we notice that the internal HAWKI configuration is actually referring to
a top-level parameter that is available to the user via the normal ``.cmds``
parameter.

If we want to use another filter, we can still use the "bang"-string format:

.. jupyter-execute::

    hawki.cmds["!OBS.filter_name"] = "H"