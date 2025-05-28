Binding ScopeSim to a local copy of the IRDB
============================================

There may be cases where we would like to have a local copy of all the files and instruments available in the Instrument Reference Database (IRDB).
In these cases it is helpful to bind out local version of ScopeSim to the IRDB.
This is quite easy and involves only two steps:

1. `Make a local copy of the IRDB`_
2. `Tell ScopeSim where to find the IRDB`_

   - `By using runtime ScopeSim commands`_
   - `By editing the ScopeSim config file`_

Make a local copy of the IRDB
-----------------------------
First we must clone the current version of the IRDB from GitHub::

    $ git clone https://github.com/AstarVienna/irdb.git

.. note:: By default this will clone the ``master`` branch of the IRDB.

To clone a specific branch, either ``checkout`` the branch you want from ``origin`` after cloning ``master`` or specify the branch at clone time with::

    $ git clone -b dev_master https://github.com/AstarVienna/irdb.git

Tell ScopeSim where to find the IRDB
------------------------------------

By using runtime ScopeSim commands
++++++++++++++++++++++++++++++++++

At the beginning of any script or jupyter notebook we can run the following command to set the default location of the instrument packages:
This is by far the simplest method, but it requires this lone of code in every script::

    import scopesim as sim
    sim.link_irdb("<path/to/IRDB>")


By editing the ScopeSim config file
++++++++++++++++++++++++++++++++++++

Find the ``defaults.yaml`` file in the top-level of the ScopeSim install directory.
This should be here::

    <path-to-python>/site-packages/scopesim/defaults.yaml

.. note:: Depending on where ScopeSim is installed, we may need administrator permissions to edit this file.

On line 36 (or there abouts), we must replace the default path string (``/inst_pkgs/``) with the path to the top level of our cloned IRDB directory (wherever we chose to put it)::

    file :
      local_packages_path : "./inst_pkgs/"  # --> change this line to point to the top level of the local IRDB folder

This path can be relative or absolute.
If a relative path is used, ScopeSim will look for instrument packages at this location relative to the directory where ScopeSim was executed (e.g. where the notebook session was started).
