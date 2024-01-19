"Generalised telescope observation simulator."

###############################################################################
#                            TURN OFF WARNINGS                                #
###############################################################################
import sys
import logging
import warnings
import yaml
from importlib import metadata
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', RuntimeWarning)    # warnings for the developer
warnings.simplefilter('ignore', category=AstropyWarning)
yaml.warnings({'YAMLLoadWarning': False})

###############################################################################
#                         PACKAGE GLOBAL VARIABLES                            #
###############################################################################

from . import rc

###############################################################################
#                         SET BASIC LOGGING LEVEL                             #
###############################################################################

# TODO: this should be replaced with YAML-based config!! see prepipy

# This should be part of ScopeSim (the app) and not scopesim_core eventually

top_logger = logging.getLogger("astar")
top_logger.setLevel(logging.WARNING)
sim_logger = top_logger.getChild(__package__)
sim_logger.setLevel(logging.DEBUG)
top_logger.propagate = False
formatter = logging.Formatter("%(name)s - %(levelname)s: %(message)s")

log_dict = rc.__config__["!SIM.logging"]
if log_dict["log_to_file"]:
    file_handler = logging.FileHandler(log_dict["file_path"],
                                       log_dict["file_open_mode"])
    file_handler.setLevel(log_dict["file_level"])  # DEBUG
    file_handler.setFormatter(formatter)
    top_logger.addHandler(file_handler)

if log_dict["log_to_console"]:
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(log_dict["console_level"])  # INFO
    stdout_handler.setFormatter(formatter)
    top_logger.addHandler(stdout_handler)


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
from .source.source import Source

from .server.database import (list_packages, download_packages, download_package,
                              list_example_data, download_example_data)

from .tests.mocks.load_basic_instrument import load_example_optical_train

###############################################################################
#                          VERSION INFORMATION                                #
###############################################################################

__version__ = metadata.version(__package__)
