Using Bang (!) strings to control ScopeSim
==========================================


TL;DR
-----

.. plot::
   :context:

   import os, scopesim
   pkg_path = os.path.join(os.getcwd(), "temp")
   scopesim.rc.__config__["!SIM.file.local_packages_path"] = pkg_path


.. plot::
   :context:
   :include-source:

   import scopesim

   cmds = scopesim.UserCommands(use_instrument="HAWKI")
   cmds["!OBS"]
   cmds["!SIM.random"]
   cmds["!SIM.random.seed"]
   cmds["!OBS.filter_name"] = "Ks"
