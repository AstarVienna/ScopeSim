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
    :class: reset, plot, format-png, format-svg

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
    format: [pdf]
    name: my_fug_A
    ---
    plt.plot([0,1], [0,1])

.. figure:: my_fug_A.png
    :name: fig:my_fug_A

    This is an included figure caption  

..
    action: [reset, clear-figure, plot]
    format: [svg]
    name: my_fug_B
    path: "./images_temp/"
    ---
    import matplotlib.pyplot as plt
    plt.plot([0,1], [1,0])

.. figure:: my_fug_B.png
    :name: fig:my_fug_B

.. code::
    :name: my_fug_C
    :class: reset, clear-figure, plot, format-png
    
    import matplotlib.pyplot as plt
    plt.plot([0,1], [1,1])

.. figure:: my_fug_C.png
    :name: fig:my_fug_C
    
    Caption to fug3

"""
