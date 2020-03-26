New features
============

Towards v0.0.3
--------------
2019-11-29
- Added ``__getitem__`` and ``__setitem__`` to both ``OpticalTrain`` and
  ``OpticsManager`` so that ``Effects`` can be accessed directly via the
  ``OpticalTrain`` by using their ``name`` value from the ``<Effect>.meta``
  dictionary. All ``Effects`` can be listed with ``<OpticalTrain>.effects``

2019-12-02
- Added support for ~/.scopesim_rc.yaml files in the user's home directory
