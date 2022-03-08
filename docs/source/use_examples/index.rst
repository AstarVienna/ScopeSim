Use-case examples
=================

A series of short and long examples illustrating ScopeSim functionality.

Users of ScopeSim are encouraged to add to this collection by either opening an
issue or submitting a pull request of GitHub.


5-liner examples
----------------

.. toctree::
    :maxdepth: 1
    :caption: Contents:
    :glob:

    5-liners/[!0]*


Longer tutorials
----------------

.. toctree::
    :maxdepth: 1
    :caption: Contents:
    :glob:

    tutorials/[!0]*



Want to add a 5-liner?
----------------------

Each use-case should ideally be a maximum of
**5 lines of ScopeSim specific code**. Obviously this not a hard and fast rule,
but simply a guideline to avoid introducing too much at once.
Import statements, setup code, etc **do not** count towards the 5 lines.

To add a 5-liner, add your code to an ipython notebook (.ipynb).
The full code goes in a ``TL;DR`` "overview" section. 
If needed, please explain what the code does in the ``Explanation`` section::

    TL;DR
    -----
    add all code here

    Explanation
    -----------
    a short explanation of what the scopesim code does

To display plots with matplotlib, make sure to include the command
``%matplotlib inline`` after ``import matplotlib``.
