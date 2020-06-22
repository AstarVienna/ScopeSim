Basic Source from Arrays
========================

TL;DR
-----

.. jupyter-execute::

    import scopesim, numpy as np

    vega = scopesim.source.source_templates.vega_spectrum(mag=20)
    point_source = scopesim.Source(x=[0], y=[0], ref=[0], spectra=[vega])
    flat_source = scopesim.Source(x=[0], y=[0], ref=[0],
                                  spectra=np.array([1, 1]),
                                  lam=np.array([0.5, 2.5]))

    print(flat_source.fields)
    print(flat_source.spectra)