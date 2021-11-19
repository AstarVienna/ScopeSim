"""The Instrument data simulator for MICADO on the E-ELT"""

################################################################################
#                            TURN OFF WARNINGS                                 #
################################################################################
import sys
import logging
import warnings
import yaml
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', RuntimeWarning)     # warnings for the developer
warnings.simplefilter('ignore', category=AstropyWarning)
yaml.warnings({'YAMLLoadWarning': False})

################################################################################
#                         PACKAGE GLOBAL VARIABLES                             #
################################################################################

from . import rc



################################################################################
#                         SET BASIC LOGGING LEVEL                              #
################################################################################

root = logging.getLogger()
root.setLevel("CRITICAL")

################################################################################
#                         IMPORT PACKAGE MODULES                               #
################################################################################

# Import all the modules to go under ScopeSim

from . import commands
from . import source
from . import optics
from . import detector
from . import effects
from . import server
from . import utils

# import specific classes from the modules to included in the global namespace

from .utils import bug_report
from .optics.optical_train import OpticalTrain
from .commands.user_commands import UserCommands
from .source.source import Source

from .server.database import list_packages, download_package

################################################################################
#                          VERSION INFORMATION                                 #
################################################################################

try:
    from .version import version as __version__
except ImportError:
    __version__ = "Version number is not available"
