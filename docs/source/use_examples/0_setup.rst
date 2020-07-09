Setup for the docs
==================

.. jupyter-execute::

    import os, scopesim

    if not os.path.exists("./temp/"):
        os.mkdir("./temp/")
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = "./temp/"
    scopesim.rc.__config__["!SIM.file.use_cached_downloads"] = False


    pkgs = ["telescopes/LFOA"] + \
           ["locations/Paranal", "telescopes/VLT", "instruments/HAWKI"] + \
           ["locations/Armazones", "telescopes/ELT", "instruments/MAORY",
            "instruments/MICADO"]

    scopesim.server.download_package(pkgs)