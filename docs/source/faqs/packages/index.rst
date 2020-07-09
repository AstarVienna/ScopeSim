Available instrument packages
=============================

.. jupyter-execute::
    :raises:
    :hide-code:
    :hide-output:

    import scopesim
    pkgs = scopesim.server.list_packages(return_pkgs=True)

To list the available packages, use ``scopesim.server.list_packages()``.

To download a package, use ``scopesim.server.download_package(<package_string>)``.
The ``package_string`` must have the following syntax: ``"<category>/<package_name>"``::

    from scopesim.server import download_package
    download_package("locations/paranal")

To download multiple packages at once, pass a list of ``package_strings``::

    download_package(["locations/Armazones",
                      "telescopes/ELT",
                      "instruments/"MICADO"])

Make sure to download all the needed packages for an instrument. E.g. for MICADO, we also need to download the ``Armazones``, ``ELT``, ``MAORY`` (optional) packages.
Running MICADO without these packages may result in erroneous simulations.

.. note:: Community contributions are most welcome!

    If you feel it could be beneficial to create an instrument package for your telescope, please feel free to do so.
    We are happy to host 3rd party packages and have them available via ScopeSim.

    To do so, please open an issue on the `IRDB GitHub page <https://github.com/astronomyk/irdb>`_ and/or submit a pull request with the new package.


Instruments
-----------
These are the currently available ``instrument`` packages

To download one or more of these packages, use the command:
``scopesim.server.download_package("instruments/<pkg_name>")``

.. jupyter-execute::
    :raises:
    :hide-code:

    for pkg in pkgs:
        if "instrument" in pkg:
            print(pkg.split("/")[1].split(".")[0])

Telescopes
----------
These are the currently available ``telescope`` packages

To download one or more of these packages, use the command:
``scopesim.server.download_package("telescopes/<pkg_name>")``

.. jupyter-execute::
    :raises:
    :hide-code:

    for pkg in pkgs:
        if "telescope" in pkg:
            print(pkg.split("/")[1].split(".")[0])

Locations
---------
These are the currently available ``location`` packages

To download one or more of these packages, use the command:
``scopesim.server.download_package("locations/<pkg_name>")``

.. jupyter-execute::
    :raises:
    :hide-code:

    for pkg in pkgs:
        if "locations" in pkg:
            print(pkg.split("/")[1].split(".")[0])


Instrument reference database
-----------------------------
All these packages are under version control on GitHub at the `Instrument Reference Database <https://github.com/astronomyk/irdb>`_. They are updated periodically.
Currently ScopeSim does not check for updates, but this functionality will be included in a later release.