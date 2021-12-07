.. image:: _static/logos/banner_red_transparent.png
    :width: 600 px
    :alt: Welcome to the ScopeSim Documentation!
    :align: center

An attempt at creating a common pythonic framework for telescope instrument
data simulators.

ScopeSim_ is on pip::

    pip install scopesim

`ScopeSim templates`_ provides templates for creating on-sky sources::

    pip install scopesim_templates

.. note:: ScopeSim only supports python 3.6 and above

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    use_examples/index
    sources/index
    effects/index
    faqs/index
    interfaces/index
    architecture/index
    Reference API <reference/modules>




Getting started
---------------
A basic simulation of ELT/MICADO image would look something like this::

    import scopesim_templates
    from scopesim.server.database import download_package
    import scopesim

    scopesim.download_package(["locations/Armazones",
                               "telescopes/ELT",
                               "instruments/MICADO"])
    cmd = scopesim.UserCommands(use_instrument="MICADO",
                                properties={"!OBS.dit": 60, "!OBS.ndit": 10,
                                            "!INST.filter_name": "Ks"})
    opt = scopesim.OpticalTrain(cmd)

    src = scopesim_templates.basic.stars.cluster()
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

Secondly we have the data describing the optical system, which are kept in
"instrument packages" held on a server. Hence the first time we want to
simulate anything, we need to use ``download_package`` to get the relevant
instrument packages from the server. In this case we want to use MICADO at the
ELT, hence we need not only the main MICADO package, but also the two support
packages for ``MICADO``: the ``ELT`` and ``Armazones``.

Finally we ``import scopesim``, the package which generates simulated
observation data sets from the two sets of input.

Once we have downloaded the packages we want to use, we need to create a
``UserCommands`` object. This contains all the information needed to
describe the model of the optical train and the how we want to observe::

    scopesim.UserCommands(use_instrument="MICADO",
                          properties={"!OBS.dit": 60, "!OBS.ndit": 10,
                                      "!INST.filter_name": "Ks"})

We start of by specifying which instrument we want to load with
``use_instrument=``. Next we pass a dictionary of keyword-value pairs
containing the settings we would like to change

Bang-strings shorten the syntax for accessing hierarchical dictionaries.
e.g. ``cmd["!OBS.ndit"]`` is the equivalent of ``cmd["OBS"]["ndit"]``


Currently available instrument packages
---------------------------------------
Below is a list of packages that are currently being maintained on our
`instrument reference database <https://github.com/astronomyk/irdb>`_.

=================== =========================== ====================================
Main Package        Support Packages            Notes
=================== =========================== ====================================
MICADO              Armazones, ELT, MAORY       Spectroscopy in Beta stage
MICADO_Sci          As above
METIS               Armazones, ELT              Only Imaging
HAWKI               Paranal, VLT
WFC3                HST                         Only NIR mode
LFOA                                            Leopold-Figl 1.5m telescope
=================== =========================== ====================================

.. warning:: We have not fully tested all packages on the latest ScopeSim

    If ScopeSim will not create an optical train with your chosen package,
    please let us know by creating a
    `Github Issue here <https://github.com/astronomyk/irdb/issues>`_.

.. note:: Community contributions welcome!

    If you think ScopeSim could be useful for your telescope/instrument,
    please don't hesitate to contact us (via email or Github) about adding a
    package to our database.


The ScopeSim python ecosystem
-----------------------------

There are several packages in the ScopeSim_ ecosystem to be aware of:

* ScopeSim_: The engine behind the whole simulator
* `ScopeSim Templates`_: A series of helper function to generate on-sky targets
* Pyckles_: Pythonic access to the Pickles (1998) spectral library and
  Brown (2014) spectral library
* IRDB_: The Instrument Reference Database, where the instrument packages are
  stored
* AnisoCADO_: For making SCAO PSF cubes that readable by ScopeSim
* skycalc_ipy_: Connects to ESOs SkyCalc server to get atmospheric spectra

.. _ScopeSim:    https://scopesim.readthedocs.io/en/latest/
.. _`ScopeSim Templates`: https://scopesim-templates.readthedocs.io/en/latest/
.. _IRDB:        https://github.com/astronomyk/irdb
.. _AnisoCADO:   https://anisocado.readthedocs.io/en/latest/
.. _skycalc_ipy: https://skycalc-ipy.readthedocs.io/en/latest/
.. _Pyckles:     https://scopesim-templates.readthedocs.io/en/latest/


.. note:: Much more information on these packages will be coming very soon!
