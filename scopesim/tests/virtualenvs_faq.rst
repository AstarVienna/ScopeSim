How to set up virtual environments for testing ScopeSim
=======================================================

1. Download all the python installers

2. Install each one in a different folder, NOT adding the python to ``PATH``
   - Make sure you also install ``Tkinter`` at the same time

3. Install ``virtualenv`` on the main python ::

    > pip install virtualenv

4. Create a virtual environment for each python. E.g::

    > cd python_virtualenvs
    > virtualenv virtualenvs\Python38 -p pythons\Python38\python.exe

5. Activate a ``virtualenv``::

    > virtualenvs\Python38\Scripts\activate

6. Install all the dependencies::

    > pip install pytest numpy scipy matplotlib astropy requests pyyaml beautifulsoup4 anisocado skycalc_ipy synphot

7. Deactivate the ``virtualenv``::

    > deactivate

8. Add ``virtualenvs`` to PyCharm settings.

- Open settings with Ctrl+Alt+S
- "Project" > "Project Interpreter"
- Click the little cog on the right side of the box "Project Interpreter"
- "Show all"
- "+"
- Add an existing interpreter by selecting the "python.exe" from the
  "Scripts" folder of your ``virtualenv`` directory. E.g::

    virtualenvs\Python38\Scripts\python.exe

