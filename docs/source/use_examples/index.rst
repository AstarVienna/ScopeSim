5-liner use-case examples
=========================

A series of quick examples illustrating ScopeSim functionality

Users or ScopeSim are encouraged to add to this collection by either opening an
issue or submitting a pull request of GitHub.

Each use-case should ideally be a maximum of
**5 lines of ScopeSim specific code**. Obviously this not a hard and fast rule,
but simply a guideline to avoid introducing too much at once.
Import statements, setup code, etc **does not** count towards the 5 lines.

.. toctree::
    :maxdepth: 1
    :caption: Contents:
    :glob:

    5-liners/*


Want to add a 5-liner?
----------------------
To add a 5-liner, add your code to a Sphinx-friendly `ReStructuredText (``.rst``)
file <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_.
The full code goes in a ``TL;DR`` "overview" section. If needed, please explain
what the code does in the ``Explanation`` section::

    TL;DR
    -----
    add all code here

    Explanation
    -----------
    a short explanation of what the scopesim code does

Code snippets are included inside a ``jupyter_execute`` tag. This code is run
in a jupyter-notebook environment and the output is also displayed::

    .. jupyter-execute::

       print("Hello World!")

.. jupyter-execute::

   print("Hello World!")

To display plots with matplotlib, make sure to include the command
``%matplotlib inline`` after ``import matplotlib``.
