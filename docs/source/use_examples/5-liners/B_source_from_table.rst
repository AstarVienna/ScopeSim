Basic Source from Astropy Tables
================================

TL;DR
-----

.. jupyter-execute::
    :raises:

    import scopesim, astropy.table as table

    vega = scopesim.source.source_templates.vega_spectrum(mag=20)
    ab_spec = scopesim.source.source_templates.ab_spectrum(mag=20)
    tbl = table.Table(names=["x",   "y",    "ref",  "weight"],
                      data=[[0, 1], [0, 2], [0, 1], [1, 0.01]])

    table_source = scopesim.Source(table=tbl, spectra=[vega, ab_spec])

    print(table_source.fields)
    print(table_source.spectra)
