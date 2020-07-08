2: Observing the same object with multiple telescopes
=====================================================
A brief introduction into using ScopeSim to observe a cluster in the LMC

.. jupyter-execute::
    :hide-code:
    :raises:

    import scopesim as sim
    sim.rc.__config__["!SIM.file.local_packages_path"] = "../temp/"


TL;DR
-----

.. jupyter-execute::

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    %matplotlib inline

    import scopesim as sim
    import scopesim_templates as sim_tp

    sim.server.download_package(["telescopes/LFOA"])
    sim.server.download_package(["locations/armazones",
                                 "telescopes/ELT",
                                 "instruments/MICADO",
                                 "instruments/MAORY"])

    cluster = sim_tp.basic.stars.cluster(mass=10000,        # Msun
                                         distance=50000,    # parsec
                                         core_radius=2,     # parsec
                                         seed=9001)         # random seed

    lfoa = sim.OpticalTrain("LFOA")
    lfoa.observe(cluster,
                 properties={"!OBS.ndit": 10, "!OBS.ndit": 360},
                 update=True)
    hdus_lfoa = lfoa.readout()

    micado = sim.OpticalTrain("MICADO")
    micado.cmds["!OBS.dit"] = 10
    micado.cmds["!OBS.ndit"] = 360
    micado.update()

    micado.observe(cluster)
    hdus_micado = micado.readout()

    plt.figure(figsize=(12,6))

    plt.subplot(121)
    plt.imshow(hdus_lfoa[0][1].data[345:385, 525:565], norm=LogNorm(), origin="lower")
    plt.colorbar()

    plt.subplot(122)
    plt.imshow(hdus_micado[0][1].data, norm=LogNorm(), origin="lower", vmax=1E6)
    plt.colorbar()
