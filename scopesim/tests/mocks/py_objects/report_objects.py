comment_plot_snippet = """
..    
    action: plot
    format: [pdf, png]
    name: my_fug
    ---
    plt.plot([0,1], [0,1])"""


literal_plot_snippet = """
.. code::
    :name: my_fug3
    :class: reset, plot, format-jpg, format-svg

    import matplotlib.pyplot as plt
    plt.plot([0,1], [1,1])"""

reset_comment_snippet = """
..
    action: reset
    ---
"""


big_rst_text = """
This parrot goes vrooom
=======================

..
    action: plot
    format: [pdf, png]
    name: my_fug
    ---
    plt.plot([0,1], [0,1])

.. figure:: my_fug.png
    :name: fig:my_fug

    This is an included figure caption  

..
    action: [reset, clear-figure, plot]
    format: [jpg, svg]
    name: my_fug2
    ---
    plt.plot([0,1], [1,0])

.. figure:: my_fug2.png
    :name: fig:my_fug2

.. code::
    :name: my_fug3
    :class: reset, clear-figure, plot, format-pdf, format-png

    plt.plot([0,1], [1,1])

.. figure:: my_fug3.png
    :name: fig:my_fug3
    
    Caption to fug3
"""
