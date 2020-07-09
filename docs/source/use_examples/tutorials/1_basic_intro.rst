1: A quick introduction to ScopeSim
===================================
A brief introduction into using ScopeSim to observe a cluster in the LMC

.. jupyter-execute::
    :hide-code:
    :raises:

    import scopesim
    scopesim.rc.__config__["!SIM.file.local_packages_path"] = "../temp/"


TL;DR
-----

.. jupyter-execute::

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    %matplotlib inline

    import scopesim as sim
    import scopesim_templates as sim_tp

    sim.server.database.download_package(["locations/Amrazones.zip",
                                          "telescopes/ELT.zip",
                                          "instruments/MAORY.zip",
                                          "instruments/MICADO.zip"])

    cluster = sim_tp.basic.stars.cluster(mass=1000,         # Msun
                                         distance=50000,    # parsec
                                         core_radius=0.3,     # parsec
                                         seed=9001)

    micado = sim.OpticalTrain("MICADO")
    micado.observe(cluster)
    hdus = micado.readout(filename="TEST.fits")

    plt.figure(figsize=(10,10))
    plt.imshow(hdus[0][1].data, norm=LogNorm(), vmax=1E5)
    plt.colorbar()
