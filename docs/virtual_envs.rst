TLDR to virtual environments
============================

1. Install the virtualenv package::

    $ pip install virtualenv

2. Install a version of Python in a folder somewhere WITHOUT adding to PATH using the normal python installers

3. Set up virtualenv with::

    $ virtualenv <name> -p <path/to/python.exe>
    $ virtualenv venv_py39 -p Python39\python.exe

4. Activate venv with::

    $ venv_py39\Scripts\activate.bat

5. Install all packages needed with pip

6. Deactivate venv with::

    $ deactivate
