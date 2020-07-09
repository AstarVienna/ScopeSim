Problem: Updating instrument packages
=====================================

Problem
-------
It doesn't appear that my packages are updating

Solution
--------
- Clear your Astropy cache:

  ``astropy.utils.data.clear_download_cache``

or

- Set the ScopeSim cache flag to False:

  ``scopesim.rc.__currsys__["!SIM.file.use_cached_downloads"] = False``

Reason
------
The packages are downloaded via the astropy ``download_file`` function and normally stored for offline use in the astropy cache.