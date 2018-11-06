Instrument Package Database
===========================

.. figure:: Inst_Pkg_Mode_Pathways.PNG
    :width: 600



.. figure:: Inst_Pkg_Server_Interfaces.PNG
    :width:


Conclusions

* The packages will be split between the instrument and telescope, but not
  between modes.
* PSF will be kept separately
* All the data for all modes for an instrument will be kept in a single
  instrument package directory
* For each of the modes there will be a config file which references the
  required files
* PSFs will be downloaded seperately from the server
* Instrument packages and telescope packages will also be downloaded separately
  as needed