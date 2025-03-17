"""Generalised telescope observation simulator."""

from importlib import metadata

###############################################################################
#                          VERSION INFORMATION                                #
###############################################################################

try:
    __version__ = metadata.version(__package__)
except metadata.PackageNotFoundError:
    __version__ = "undetermined"


###############################################################################
#                            TURN OFF WARNINGS                                #
###############################################################################

import warnings
import yaml
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', RuntimeWarning)  # warnings for the developer
warnings.simplefilter('default', DeprecationWarning)  # allow in general
warnings.simplefilter('ignore', category=AstropyWarning)
yaml.warnings({'YAMLLoadWarning': False})

###############################################################################
#                         PACKAGE GLOBAL VARIABLES                            #
###############################################################################

from . import rc

###############################################################################
#                         SET BASIC LOGGING LEVEL                             #
###############################################################################

# This should be part of ScopeSim (the app) and not scopesim_core eventually
# TODO: need to add a function to reload the config!

# Import convenience functions
from .utils import update_logging, log_to_file, set_console_log_level
update_logging()


###############################################################################
#                         IMPORT PACKAGE MODULES                              #
###############################################################################

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
from .commands.scopesimple import Simulation
from .source.source import Source

from .server.database import (list_packages, download_packages, download_package,
                              list_example_data, download_example_data)

from .tests.mocks.load_basic_instrument import load_example_optical_train, example_simulation
