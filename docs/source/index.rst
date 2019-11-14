Welcome to ScopeSim's documentation!
====================================

An attempt at creating a common pythonic framework for telescope instrument
data  simulators.

ScopeSim_ hasn't yet been released on pip, but it is in a usable state::

    pip install git+https://github.com/astronomyk/ScopeSim

The same applies for the templates package: `ScopeSim templates`_::

    pip install git+https://github.com/astronomyk/ScopeSim_templates

.. note:: ScopeSim only supports python 3.5 and above

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    interfaces/controlling_simulations
    examples/index
    interfaces/index
    architecture/index
    effects/index


Getting started
---------------
A basic simulation of VLT/HAWKI image would look something like this::

    import scopesim_templates
    from scopesim.server.database import download_package
    import scopesim

    download_package(["Paranal", "VLT", "HAWKI"])
    cmd = scopesim.UserCommands(use_instrument="HAWKI",
                                properties={"!OBS.dit": 60, "!OBS.ndit": 10,
                                            "!INST.filter_name": "Ks"})
    opt = scopesim.OpticalTrain(cmd)

    src = scopesim_templates.stars.open_cluster()
    opt.observe(src)
    opt.readout().writeto("my_image.fits")

Let's break this down a bit.

There are three major components of any simulation workflow:

1. the target description,
2. the telescope/instrument model, and
3. the observation.

The three ``import`` statements are representative of this.
Firstly we have ``import scopesim_templates`` which the package that can
generate a description of our target. It contains lots of helper functions
which return the ``Source`` objects accepted during a ``scopesim`` observation.
Secondly we have the data descriping the optical system, which are kept in
"instrument packages" held on a server. Hence the first time we want to
simulate anything, we need to use ``download_package`` to get the relevant
instrument packages from the server. In this case we want to use HAWKI at the
VLT, hence we need not only the main HAWKI package, but also the two support
packages for ``HAWKI``: the ``VLT`` and ``Paranal``.
Finally we ``import scopesim``, the package which generates simulated
observation data sets from the two sets of input.

Once we have downloaded the packages we want to use, we need to create a
``UserCommands`` object. This contains all the information needed to
describe the model of the optical train and the how we want to observe::

    scopesim.UserCommands(use_instrument="HAWKI",
                          properties={"!OBS.dit": 60, "!OBS.ndit": 10,
                                      "!INST.filter_name": "Ks"})

We start of by specifying which instrument we want to load with
``use_instrument=``. Next we pass a dictionary of keyword-value pairs
containing the settings we would like to change


Bang-strings


The ScopeSim python ecosystem
-----------------------------

There are several packages in the ScopeSim_ ecosystem to be aware of:

* ScopeSim_: The engine behind the whole simulator
* `ScopeSim Templates`_: A series of helper function to generate on-sky targets
* IRDB_: The Instrument Reference Database, where the instrument packages are
  stored
* AnisoCADO_: For making SCAO PSF cubes that readable by ScopeSim
* skycalc_ipy_: Connects to ESOs SkyCalc server to get atmospheric spectra

.. _ScopeSim: https://github.com/astronomyk/ScopeSim
.. _`ScopeSim Templates`: https://github.com/astronomyk/ScopeSim_Templates
.. _IRDB: https://github.com/astronomyk/irdb
.. _AnisoCADO: https://github.com/astronomyk/AnisoCADO
.. _skycalc_ipy: https://github.com/astronomyk/skycalc_ipy




.. note:: Much more information on these packages will be coming very soon!





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
