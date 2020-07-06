Using Bang (!) strings to control ScopeSim
==========================================

TL;DR
-----

.. jupyter-execute::
    :hide-code:
    :raises:

    import os, scopesim
    pkg_path = os.path.join(os.getcwd(), "temp")
    if not os.path.exists(pkg_path):
        os.mkdir(pkg_path)
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = pkg_path

.. jupyter-execute::
    :raises:

    hawki = scopesim.OpticalTrain("HAWKI")

    hawki.effects
    hawki["detector_linearity"].include = False
    hawki["detector_linearity"].meta["include"] = True


Background
----------

This 5-liner uses concepts from:
- :doc:`loading_packages`: downloading instrument packages
- :doc:`bang_strings`: accessing top- and lower-level parameters.


Explanation
-----------

To list all the effects in the HAWKI optical train, we do:

.. jupyter-execute::

    hawki = scopesim.OpticalTrain("HAWKI")
    hawki.effects

This table already shows us which ``Effect`` objects are turned on.

To turn ``Effect`` object on or off manually, we use the ``.include`` attribute.
Here we must call the ``Effect`` by it's name as given in the previous table:

.. jupyter-execute::

    hawki["detector_linearity"].include = False
    hawki["detector_linearity"].include

Turning back on is simple:

.. jupyter-execute::

    hawki["detector_linearity"].include = True

If we want to change many parameters at once, including ``include`` we can
access it via the ``.meta`` dictionary:

.. jupyter-execute::

    hawki["detector_linearity"].meta["include"] = True