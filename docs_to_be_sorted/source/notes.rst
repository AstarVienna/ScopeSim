* Making effects from YAMLs
    * The ``properties`` dictionary from any high level OpticalElement will be
      added to the ``kwargs`` of any Effect as it is created

* Too slow due to ridiculous amounts of FOVs
    * The wavelength range is too large. Make sure the the Filter effect is
      passing ``minimum_throughput`` larger than the wing transmission. I.e. if
      the TC has 0.0001 wings, set ``minimum_throughput : 0.00011`` in the
      Filter effect ``kwargs`` section
    * Increase the chunk size by setting ``!SIM.computing.chunk_size``
