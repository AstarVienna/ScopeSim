How to add extensions to Sphinx
===============================

1. Add the name of the extension module (i.e. ``scopesim_sphinx_ext``) to ``conf.py`` in the ``extensions`` list
2. Add the path of the extensions folder (i.e. ``_ext``) to the python path at the beginning of ``conf.py``:

    sys.path.append(os.path.abspath("./_ext"))

3. Add a setup function to the extension module (i.e. ``scopesim_sphinx_ext``):

    def setup(app):
        with open("test.txt", "w") as f:
            f.write("Hello World!")

4. The base directory for running this is where ``make html`` is run (i.e. ``docs/``)
