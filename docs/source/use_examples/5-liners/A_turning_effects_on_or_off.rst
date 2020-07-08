Control: Turning Effect objects on or off
=========================================

.. jupyter-execute::
    :hide-code:
    :raises:

    import scopesim
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = "./temp/"

TL;DR
-----

.. jupyter-execute::
    :raises:

    hawki = scopesim.OpticalTrain("HAWKI")

    hawki.effects
    hawki["detector_linearity"].include = False
    hawki["detector_linearity"].meta["include"] = True


Background
----------

This 5-liner uses concepts from:

- :doc:`A_loading_packages`: downloading instrument packages
- :doc:`A_bang_strings`: accessing top- and lower-level parameters.


Explanation
-----------

To list all the effects in the HAWKI optical train, we do:

.. jupyter-execute::
    :raises:

    hawki = scopesim.OpticalTrain("HAWKI")
    print(hawki.effects)

This table already shows us which ``Effect`` objects are turned on.

To turn ``Effect`` object on or off manually, we use the ``.include`` attribute.
Here we must call the ``Effect`` by it's name as given in the previous table:

.. jupyter-execute::
    :raises:

    hawki["detector_linearity"].include = False
    hawki["detector_linearity"].include

Turning back on is simple:

.. jupyter-execute::
    :raises:

    hawki["detector_linearity"].include = True

If we want to change many parameters at once, including ``include`` we can
access it via the ``.meta`` dictionary:

.. jupyter-execute::
    :raises:

    hawki["detector_linearity"].meta["include"] = True
