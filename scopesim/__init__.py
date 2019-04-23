"""The Instrument data simulator for MICADO on the E-ELT"""

################################################################################
#                            TURN OFF WARNINGS                                 #
################################################################################

import warnings
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', RuntimeWarning)     # warnings for the developer
warnings.simplefilter('ignore', category=AstropyWarning)


################################################################################
#                         PACKAGE GLOBAL VARIABLES                             #
################################################################################

from . import rc

################################################################################
#                         ENABLE/DISABLE LOGGING                               #
################################################################################

import logging

if rc.__rc__["SIM_LOGGING"]:
    logging.basicConfig(filename=rc.__rc__["SIM_LOGGING_FILE"],
                        filemode='w',
                        level=logging.getLevelName(
                            rc.__rc__["SIM_LOGGING_LEVEL"]),
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info("ScopeSim imported, logging initiated")


################################################################################
#                         IMPORT PACKAGE MODULES                               #
################################################################################

# Import all the modules to go under ScopeSim
from . import utils

from . import spectral
from . import spatial
from scopesim.optics import OLD_psf

from . import detector
from . import optics
from . import commands
from .source import source

from . import simulation

# import specific Classes from the modules to be accessible in the global
# namespace
from .utils import bug_report
from .utils import get_extras


################################################################################
#                          VERSION INFORMATION                                 #
################################################################################

try:
    from .version import version as __version__
except ImportError:
    __version__ = "Version number is not available"
