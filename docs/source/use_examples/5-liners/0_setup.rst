Setup for the docs
==================

.. jupyter-execute::

    import os, scopesim

    if not os.path.exists("./temp/"):
        os.mkdir("./temp/")
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = "./temp/"
    pkg_names = ["locations/Paranal", "telescopes/VLT", "instruments/HAWKI"]
    scopesim.server.download_package(pkg_names)
