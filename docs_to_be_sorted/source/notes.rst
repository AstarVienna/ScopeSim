* Making effects from YAMLs
    * The ``properties`` dictionary from any high level OpticalElement will be
      added to the ``kwargs`` of any Effect as it is created

* Too slow due to ridiculous amounts of FOVs
    * The wavelength range is too large. Make sure the the Filter effect is
      passing ``minimum_throughput`` larger than the wing transmission. I.e. if
      the TC has 0.0001 wings, set ``minimum_throughput : 0.00011`` in the
      Filter effect ``kwargs`` section
    * Increase the chunk size by setting ``!SIM.computing.chunk_size``

* Access effects directly with names
    * When an OpticalTrain is created, the Effect objects are also initialised
      If you want to change values directly in an Effect object, you can access
      them by their names::

        >>> opt = scopesim.OpticalTrain("MICADO")
        >>> opt.OpticsManager
        OpticsManager contains 7 OpticalElements
        [0] "misc" contains 0 effects
        [1] "armazones" contains 2 effects
        [2] "ELT" contains 2 effects
        ...
        >>> opt.optics_manager["MICADO"]
        OpticalElement : "MICADO" contains 3 Effects:
        [0] SurfaceList: "micado_static_surfaces"
        [1] FilterCurve: "micado_filter"
        ...
        >>> opt.optics_manager["MICADO"]["micado_filter"]
        FilterCurve: "micado_filter"

      The data relevant to the flux through the OpticalTrain are kept in the
      ``.data`` attribute::

        >>> opt.optics_manager["MICADO"]["micado_filter"].data
        <Table length=1718>
        wavelength transmission
         float64     float64
        ---------- ------------
             0.783       0.0001
             0.784       0.0001
        ...

      The useful data of an Effect object is contained in the ``.meta``
      dictionary attribute::

        >>> opt.optics_manager["MICADO"]["micado_filter"].meta
        ...
        'center': 2.144876984095012,
        'width': 0.3470169216021642,
        ...

      .. attention::
         Often the settings for how to use the data will also be stored in the
         ``.meta`` attribute. This is a normal dictionary, and changing values
         here can cause different (possibly un-expected) behaviour at run-time.
         Use this functionality with care.

         The recommended way is to change the info in the ``UserCommands``
         object which spawned the ``OpticalTrain`` object, then reload the
         ``OpticalTrain``.