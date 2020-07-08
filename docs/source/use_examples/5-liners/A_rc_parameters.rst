Control: Global simulation parameters
=====================================

TL;DR
-----

.. jupyter-execute::
    :raises:

    import scopesim

    scopesim.rc.__currsys__["!SIM.random.seed"] = 9001
    scopesim.rc.__currsys__["!SIM.file.local_packages_path"] = "./"

    print(scopesim.rc.__currsys__["!SIM"])
    print(scopesim.rc.__currsys__["!SIM.file"])
