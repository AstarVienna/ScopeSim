## Author : Kieran
## Date-created :  2018-11-09
## Date-modified : 2018-11-09
## Changes : 
##    2018-11-09 Added keywords INST_PKG_LOCAL_PATH and INST_PKG_SERVER_PATH 
##
## This file contains the default parameters regarding the software,
## It does NOT contain information regarding the instrument or telescope 
## configuration

import os
import inspect

# Package directory
PKG_DIR = os.path.dirname(inspect.getfile(inspect.currentframe()))

# Instrument package paths
INST_PKG_LOCAL_PATH     = "inst_pkgs/"
INST_PKG_LOCAL_DB_NAME  = "LocalInstPkgDB.dat"
INST_PKG_SERVER_PATH    = "https://www.univie.ac.at/simcado/InstPkgSvr/"
INST_PKG_SERVER_DB_NAME = "InstPkgDB.dat"
